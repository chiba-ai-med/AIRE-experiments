from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("8.10.0")

configfile: "workflow/config.yaml"

DOCKER_IMAGE = "koki/aire_workflow:dev"
CONDA_NF = "envs/nextflow.yaml"

#################################
# SRA Hi-C: SRA fastq --> nf-core/hic --> bulk mTEC mcool
#################################
# Why this is a separate Snakefile from preprocess.smk:
# - Heavy compute (~30 GB SRA download, ~200 GB intermediate, 12-24 h on 32 cores)
# - Wraps two nf-core Nextflow pipelines (fetchngs + hic) which manage their
#   own conda/docker envs internally, so we only need nextflow + java in the
#   surrounding env.
# - Failure modes (network, SRA throttling, nf-core version drift) are
#   different from the main pipeline; isolate them.
#
# Pipeline:
#   write_sra_ids       : config sra_samples -> sra_ids.csv (one SRX per line)
#   fetchngs            : nf-core/fetchngs SRX list -> fastq.gz + samplesheet.csv
#   run_nfcore_hic      : nf-core/hic samplesheet -> per-sample mcool
#   merge_bulk_mcool    : sum config bulk_replicates -> data/mtec/processed/mtec_hic.mcool
#
# The final mtec_hic.mcool feeds back into preprocess.smk's bin_hic_by_chr
# rule (wildcard tissue=mtec), which then becomes Machima2 input.
#
# Usage:
#   With conda:        bash workflow/run_sra_hic.sh
#   With Singularity:  bash workflow/run_sra_hic.sh --use-singularity

rule all:
    input:
        'data/mtec/processed/mtec_hic.mcool'

#################################
# Write SRX accession list for nf-core/fetchngs
#################################
rule write_sra_ids:
    """Emit one SRX per line; consumed by nf-core/fetchngs as --input."""
    output:
        ids='data/mtec/sra/sra_ids.csv'
    run:
        from pathlib import Path
        Path(output.ids).parent.mkdir(parents=True, exist_ok=True)
        with open(output.ids, 'w') as f:
            for sample_name, srx in config['mtec']['query']['sra_samples'].items():
                f.write(f"{srx}\n")

#################################
# nf-core/fetchngs: SRA -> fastq.gz + samplesheet
#################################
rule fetchngs:
    """
    Run nf-core/fetchngs to resolve each SRX into its SRRs and download
    paired fastq.gz. Produces a samplesheet directly compatible with
    nf-core/hic. Internally handles per-SRR prefetch + fasterq-dump.
    """
    input:
        ids='data/mtec/sra/sra_ids.csv'
    output:
        samplesheet='data/mtec/sra/fetchngs/samplesheet/samplesheet.csv'
    params:
        outdir='data/mtec/sra/fetchngs',
        # fetchngs maps download_method options: aspera | sratools | ftp.
        # sratools is the most reliable default; aspera is faster if available.
        download_method='sratools'
    log:
        'logs/sra_hic/fetchngs.log'
    benchmark:
        'benchmarks/sra_hic/fetchngs.txt'
    conda:
        CONDA_NF
    container:
        f"docker://{DOCKER_IMAGE}"
    threads: 8
    shell:
        # NOTE: --nf_core_pipeline is intentionally omitted. As of fetchngs
        # 1.12.0 the supported choices are rnaseq | atacseq | viralrecon |
        # taxprofiler -- "hic" is not one of them. Without the flag, fetchngs
        # writes a generic samplesheet whose sample/fastq_1/fastq_2 columns
        # are directly compatible with nf-core/hic's input schema.
        'nextflow run nf-core/fetchngs '
        '-r 1.12.0 '
        '-profile conda '
        '--input {input.ids} '
        '--outdir {params.outdir} '
        '--download_method {params.download_method} '
        '-resume '
        '>& {log}'

#################################
# Remap fetchngs samplesheet sample column (SRX -> friendly name)
#################################
rule remap_samplesheet:
    """
    fetchngs writes the samplesheet with sample = SRX accession because
    --nf_core_pipeline supports only rnaseq/atacseq/viralrecon/taxprofiler
    (not hic). We rewrite the sample column using the SRX -> friendly-name
    mapping in config.mtec.query.sra_samples so that downstream nf-core/hic
    output cool files use the friendly names referenced in bulk_replicates.
    """
    input:
        samplesheet='data/mtec/sra/fetchngs/samplesheet/samplesheet.csv'
    output:
        samplesheet='data/mtec/sra/fetchngs/samplesheet/samplesheet_renamed.csv'
    run:
        import csv
        srx_to_name = {srx: name for name, srx in config['mtec']['query']['sra_samples'].items()}
        with open(input.samplesheet) as f, open(output.samplesheet, 'w', newline='') as g:
            reader = csv.reader(f)
            writer = csv.writer(g, quoting=csv.QUOTE_MINIMAL)
            header = next(reader); writer.writerow(header)
            sidx = header.index('sample')
            for row in reader:
                row[sidx] = srx_to_name.get(row[sidx], row[sidx])
                writer.writerow(row)

#################################
# nf-core/hic: fastq -> per-sample mcool
#################################
rule run_nfcore_hic:
    """
    Run nf-core/hic with the renamed samplesheet. Outputs per-sample
    multi-resolution cools at the resolutions in config.resolutions
    plus the project base 1000 bp.

    Parameters are passed via -params-file (YAML) rather than CLI flags
    because Nextflow 25.x's nf-validation strictly types CLI values --
    --bin_size 100000 is parsed as Integer and rejected against the
    schema's String type. YAML preserves the quoted string type.
    """
    input:
        samplesheet='data/mtec/sra/fetchngs/samplesheet/samplesheet_renamed.csv'
    output:
        sentinel='data/mtec/processed/nfcore_hic/.done',
        params_yaml='data/mtec/processed/nfcore_hic/params.yaml'
    params:
        outdir='data/mtec/processed/nfcore_hic',
        bin_sizes=','.join(str(r) for r in config['resolutions']),
        restriction_site=config['mtec']['query']['restriction_site'],
        ligation_site=config['mtec']['query']['ligation_site'],
        genome=config['genome']
    log:
        'logs/sra_hic/nfcore_hic.log'
    benchmark:
        'benchmarks/sra_hic/nfcore_hic.txt'
    conda:
        CONDA_NF
    container:
        f"docker://{DOCKER_IMAGE}"
    threads: 16
    shell:
        'mkdir -p {params.outdir} && '
        'printf "%s\\n" '
        '  \'input: "{input.samplesheet}"\' '
        '  \'outdir: "{params.outdir}"\' '
        '  \'genome: "{params.genome}"\' '
        '  \'bin_size: "{params.bin_sizes}"\' '
        '  \'restriction_site: "{params.restriction_site}"\' '
        '  \'ligation_site: "{params.ligation_site}"\' '
        '  \'skip_multiqc: true\' '
        '  > {output.params_yaml} && '
        'nextflow run nf-core/hic '
        '-r 2.1.0 '
        '-profile conda '
        '-c workflow/nfcore_hic.config '
        '-params-file {output.params_yaml} '
        '-resume '
        '>& {log} '
        '&& touch {output.sentinel}'

#################################
# Merge WT replicates into the bulk mTEC mcool (Machima2 X_Epi)
#################################
rule merge_bulk_mcool:
    """
    Sum the per-sample cools listed in config.bulk_replicates at the
    finest resolution, then zoomify+balance to the full multi-resolution
    mcool used as Machima2 X_Epi for mTEC. Per-sample (KO included) cools
    remain available under data/mtec/processed/nfcore_hic/ for evaluation.
    """
    input:
        sentinel='data/mtec/processed/nfcore_hic/.done'
    output:
        mcool='data/mtec/processed/mtec_hic.mcool'
    params:
        nfcore_outdir='data/mtec/processed/nfcore_hic',
        replicates=lambda wc: config['mtec']['query']['bulk_replicates'],
        base_resolution=min(config['resolutions']),
        zoom_resolutions=','.join(str(r) for r in sorted(set(config['resolutions']),
                                                          reverse=True))
    log:
        'logs/sra_hic/merge_bulk_mcool.log'
    benchmark:
        'benchmarks/sra_hic/merge_bulk_mcool.txt'
    conda:
        CONDA_NF
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'bash src/sra_merge_bulk_mcool.sh '
        '"{params.nfcore_outdir}" "{params.replicates}" '
        '{params.base_resolution} "{params.zoom_resolutions}" {output.mcool} '
        '>& {log}'
