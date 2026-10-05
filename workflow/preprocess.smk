# Preprocess phase: data -> per-chromosome bin matrices.
#
# Brain (practice):
#   - reference: GSE210747 (10x multiome scRNA+scATAC, mouse cortex E15.5)
#       fragments + DGEM TSV --> MACS2 peaks --> .h5mu --> leiden labels
#   - query    : GSE96107 (bulk Hi-C, FACS-purified ncx NPC + cortical neurons)
#       misha tracks --> bg2 --> mcool (combined + per-celltype) --> bin matrices
# mTEC (production):
#   - reference: local cellranger-arc output at data/mtec/reference
#       atac_fragments.tsv.gz + filtered matrix --> MACS2 peaks --> .h5mu
#   - query    : NOT in this Snakefile -- see workflow/sra_hic.smk
#
# All shared constants (CONDA_*, DOCKER_IMAGE, CHROM_SIZES, MISHA_GDB,
# DEFAULT_RESOLUTION) live in common.smk and are inherited at include time.

#################################
# Phase target rule.
#################################
rule preprocess:
    """All preprocess outputs needed before run_machima2 can fire."""
    input:
        # GEO downloads
        'data/brain/reference/GSE210747/.downloaded',
        'data/brain/query/GSE96107/.downloaded',
        # MACS2 peaks
        'data/brain/processed/brain_peaks.narrowPeak',
        'data/mtec/processed/mtec_peaks.narrowPeak',
        # Multiome MuData (raw + leiden-labelled)
        'data/brain/processed/multiome.h5mu',
        'data/mtec/processed/multiome.h5mu',
        'data/brain/processed/multiome_labeled.h5mu',
        'data/mtec/processed/multiome_labeled.h5mu',
        # Hi-C bulk (combined) + per-celltype mcools
        'data/brain/processed/brain_hic.mcool',
        'data/brain/processed/brain_hic_npc.mcool',
        'data/brain/processed/brain_hic_cn.mcool',
        # Per-chromosome bins (ATAC + Hi-C combined + Hi-C per-celltype)
        f'data/brain/processed/atac_bins_{DEFAULT_RESOLUTION}/',
        f'data/mtec/processed/atac_bins_{DEFAULT_RESOLUTION}/',
        f'data/brain/processed/hic_bins_{DEFAULT_RESOLUTION}/',
        f'data/brain/processed/hic_bins_npc_{DEFAULT_RESOLUTION}/',
        f'data/brain/processed/hic_bins_cn_{DEFAULT_RESOLUTION}/',
        f'data/mtec/processed/hic_bins_{DEFAULT_RESOLUTION}/',
        # leiden cluster -> celltype map (consumed by validate_against_sorted_bulk)
        'data/brain/processed/cluster_celltype_map.tsv',

#################################
# Download: brain reference (10x multiome supplementary)
#################################
rule download_brain_reference:
    """
    Download GSE210747 series-level supplementary files into
    data/brain/reference/GSE210747/. Includes per-condition
    atac_fragments.tsv.gz and the combined RNA DGEM TSV.
    """
    output:
        sentinel='data/brain/reference/GSE210747/.downloaded'
    params:
        accession='GSE210747',
        outdir='data/brain/reference/GSE210747'
    log:
        'logs/preprocess/download_brain_reference.log'
    benchmark:
        'benchmarks/preprocess/download_brain_reference.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'bash src/preprocess_download_geo.sh series {params.accession} {params.outdir} >& {log} '
        '&& touch {output.sentinel}'

#################################
# Download: brain query (Bonev Hi-C ncx samples)
#################################
rule download_brain_query:
    """
    Download per-GSM Hi-C pairs tarballs for the FACS-purified ncx
    samples (NPC + cortical neurons, ~141 GB total). Stored under
    data/brain/query/GSE96107/<GSM>/.
    """
    output:
        sentinel='data/brain/query/GSE96107/.downloaded'
    params:
        gsm_list=lambda wc: ' '.join(config['brain']['query']['samples']),
        outdir='data/brain/query/GSE96107'
    log:
        'logs/preprocess/download_brain_query.log'
    benchmark:
        'benchmarks/preprocess/download_brain_query.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'bash src/preprocess_download_geo.sh samples "{params.gsm_list}" {params.outdir} >& {log} '
        '&& touch {output.sentinel}'

#################################
# Peak calling: brain primary-sample ATAC fragments
#################################
rule call_peaks_brain:
    """
    MACS2 narrowPeak calling on the brain primary-sample fragments
    (config.brain.reference.primary_sample). Standard ATAC --shift -75
    --extsize 150 recipe; cell-barcode column is ignored at peak-call time
    and re-attached when the cell x peak matrix is built.
    """
    input:
        sentinel='data/brain/reference/GSE210747/.downloaded'
    output:
        peaks='data/brain/processed/brain_peaks.narrowPeak'
    params:
        fragments=f"data/brain/reference/{config['brain']['reference']['accession']}/"
                  f"{config['brain']['reference']['accession']}_"
                  f"{config['brain']['reference']['primary_sample']}_atac_fragments.tsv.gz",
        outdir='data/brain/processed',
        name='brain'
    log:
        'logs/preprocess/call_peaks_brain.log'
    benchmark:
        'benchmarks/preprocess/call_peaks_brain.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'bash src/preprocess_call_peaks.sh {params.fragments} {params.name} {params.outdir} >& {log}'

#################################
# Peak calling: mTEC ATAC fragments
#################################
rule call_peaks_mtec:
    """
    MACS2 narrowPeak calling on cellranger-arc atac_fragments.tsv.gz.
    cellranger's own peaks (atac_peaks.bed) are intentionally not used --
    we re-call peaks here so brain and mTEC share the exact same peak-set
    construction (see project_multiome_input_unification memory).
    """
    input:
        fragments='data/mtec/reference/atac_fragments.tsv.gz'
    output:
        peaks='data/mtec/processed/mtec_peaks.narrowPeak'
    params:
        outdir='data/mtec/processed',
        name='mtec'
    log:
        'logs/preprocess/call_peaks_mtec.log'
    benchmark:
        'benchmarks/preprocess/call_peaks_mtec.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'bash src/preprocess_call_peaks.sh {input.fragments} {params.name} {params.outdir} >& {log}'

#################################
# Build multiome MuData: brain
#################################
rule build_mudata_brain:
    """
    Compose data/brain/processed/multiome.h5mu from the combined RNA DGEM
    TSV and the primary-sample ATAC fragments + MACS2 peaks. Cell barcodes
    are intersected, which auto-subsets the combined DGEM to the cells
    from the primary protocol condition.
    """
    input:
        sentinel='data/brain/reference/GSE210747/.downloaded',
        peaks='data/brain/processed/brain_peaks.narrowPeak'
    output:
        h5mu='data/brain/processed/multiome.h5mu'
    params:
        rna=f"data/brain/reference/{config['brain']['reference']['accession']}/"
            f"{config['brain']['reference']['accession']}_Mouse_cortex_DGEM.tsv.gz",
        atac_fragments=f"data/brain/reference/{config['brain']['reference']['accession']}/"
                       f"{config['brain']['reference']['accession']}_"
                       f"{config['brain']['reference']['primary_sample']}_atac_fragments.tsv.gz",
        primary_sample=config['brain']['reference']['primary_sample'],
        genome=config['genome']
    log:
        'logs/preprocess/build_mudata_brain.log'
    benchmark:
        'benchmarks/preprocess/build_mudata_brain.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/preprocess_build_multiome_mudata.py '
        '--rna-format dgem-tsv '
        '--rna-input {params.rna} '
        '--primary-sample {params.primary_sample} '
        '--atac-fragments {params.atac_fragments} '
        '--atac-peaks {input.peaks} '
        '--genome {params.genome} '
        '--output {output.h5mu} >& {log}'

#################################
# Build multiome MuData: mTEC
#################################
rule build_mudata_mtec:
    """
    Compose data/mtec/processed/multiome.h5mu from the cellranger-arc 10x
    mtx (Gene Expression feature type only) + atac_fragments.tsv.gz +
    MACS2 peaks. cellranger's own peaks/peak-counts are not used.
    """
    input:
        rna_mtx='data/mtec/reference/filtered_feature_bc_matrix/matrix.mtx.gz',
        atac_fragments='data/mtec/reference/atac_fragments.tsv.gz',
        peaks='data/mtec/processed/mtec_peaks.narrowPeak'
    output:
        h5mu='data/mtec/processed/multiome.h5mu'
    params:
        rna_dir='data/mtec/reference/filtered_feature_bc_matrix',
        genome=config['genome']
    log:
        'logs/preprocess/build_mudata_mtec.log'
    benchmark:
        'benchmarks/preprocess/build_mudata_mtec.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/preprocess_build_multiome_mudata.py '
        '--rna-format 10x-mtx '
        '--rna-input {params.rna_dir} '
        '--atac-fragments {input.atac_fragments} '
        '--atac-peaks {input.peaks} '
        '--genome {params.genome} '
        '--output {output.h5mu} >& {log}'

#################################
# Cluster + label: brain
#################################
rule cluster_label_brain:
    """
    Leiden cluster the RNA modality and write per-cell labels into
    mdata.obs['cell_type'] as 'cluster_<k>' strings. For multiome the
    barcode is shared so the label transfers to ATAC by index alignment
    (no explicit transfer step needed).
    """
    input:
        h5mu='data/brain/processed/multiome.h5mu'
    output:
        h5mu='data/brain/processed/multiome_labeled.h5mu'
    log:
        'logs/preprocess/cluster_label_brain.log'
    benchmark:
        'benchmarks/preprocess/cluster_label_brain.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/preprocess_cluster_label.py '
        '--input {input.h5mu} --output {output.h5mu} >& {log}'

#################################
# Cluster + label: mTEC
#################################
rule cluster_label_mtec:
    """As cluster_label_brain, on the mTEC multiome."""
    input:
        h5mu='data/mtec/processed/multiome.h5mu'
    output:
        h5mu='data/mtec/processed/multiome_labeled.h5mu'
    log:
        'logs/preprocess/cluster_label_mtec.log'
    benchmark:
        'benchmarks/preprocess/cluster_label_mtec.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/preprocess_cluster_label.py '
        '--input {input.h5mu} --output {output.h5mu} >& {log}'

#################################
# Shared: write <genome>.chrom.sizes (single source for misha + cooler + binning)
#################################
rule write_chrom_sizes:
    """
    Dump <genome>.chrom.sizes from snapatac2.genome.<genome>. The same file
    feeds the misha gdb setup (extract_bonev_track), cooler load/zoomify
    (bg2_to_mcool), and any chr-restricted downstream step.
    """
    output:
        sizes=CHROM_SIZES
    params:
        genome=config['genome']
    log:
        'logs/preprocess/write_chrom_sizes.log'
    benchmark:
        'benchmarks/preprocess/write_chrom_sizes.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/preprocess_write_chrom_sizes.py '
        '--genome {params.genome} --output {output.sizes} >& {log}'

#################################
# One-shot: download misha genome database (~786 MB)
#################################
rule setup_misha_genome:
    """
    Download the misha pre-built genome database via gdb.create_genome().
    Persistent cache: <MISHA_GDB>/. Subsequent rules (extract_bonev_track)
    just call gsetroot() against it -- no re-download.
    """
    output:
        sentinel=f"{MISHA_GDB}/.created"
    params:
        genome=config['genome'],
        parent='data/misha_genomes'
    log:
        'logs/preprocess/setup_misha_genome.log'
    benchmark:
        'benchmarks/preprocess/setup_misha_genome.txt'
    conda:
        CONDA_R
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'Rscript src/preprocess_setup_misha_genome.R '
        '{params.genome} {params.parent} >& {log}'

#################################
# Bonev misha track -> BedGraph2D (per GSM)
#################################
rule extract_bonev_track_brain:
    """
    Extract a Bonev/Tanay misha 2D track tarball into a BedGraph2D file
    at the project's finest bin resolution (min(config.resolutions)). The
    Bonev tarballs ship StatQuadTreeCached binary 2D tracks that only the
    misha R package can deserialise; this rule wraps that and emits a
    standard bg2.gz that cooler load can ingest.

    Tarball filename pattern is `<GSM>_<sample>.tar.gz` inside
    `data/brain/query/GSE96107/<GSM>/`; the inner `<sample>` is resolved
    via glob at runtime.
    """
    input:
        sentinel='data/brain/query/GSE96107/.downloaded',
        misha_ready=f"{MISHA_GDB}/.created"
    output:
        bg2='data/brain/processed/bonev_bg2/{gsm}.bg2.gz'
    wildcard_constraints:
        gsm=r'GSM\d+'
    params:
        tarball=lambda wc: glob.glob(f'data/brain/query/GSE96107/{wc.gsm}/{wc.gsm}_*.tar.gz')[0],
        misha_root=MISHA_GDB,
        resolution=min(config['resolutions'])
    log:
        'logs/preprocess/extract_bonev_track_brain_{gsm}.log'
    benchmark:
        'benchmarks/preprocess/extract_bonev_track_brain_{gsm}.txt'
    conda:
        CONDA_R
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'Rscript src/preprocess_extract_bonev_track.R '
        '{params.tarball} {params.misha_root} {params.resolution} {output.bg2} '
        '>& {log}'

#################################
# Aggregate per-GSM bg2 into bulk multi-resolution mcool (brain, combined)
#################################
rule bg2_to_mcool_brain:
    """
    Aggregate per-GSM BedGraph2D files (NPC + CN reps) into a single
    bulk multi-resolution mcool. cooler load -f bg2 per-GSM, then
    cooler merge sums counts at coincident bins, then zoomify+balance.
    """
    input:
        bg2_files=expand('data/brain/processed/bonev_bg2/{gsm}.bg2.gz',
                         gsm=config['brain']['query']['samples']),
        chrom_sizes=CHROM_SIZES
    output:
        mcool='data/brain/processed/brain_hic.mcool'
    params:
        base_resolution=min(config['resolutions']),
        zoom_resolutions=','.join(str(r) for r in sorted(set(config['resolutions']),
                                                          reverse=True))
    log:
        'logs/preprocess/bg2_to_mcool_brain.log'
    benchmark:
        'benchmarks/preprocess/bg2_to_mcool_brain.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'bash src/preprocess_bg2_to_mcool.sh '
        '{input.chrom_sizes} {params.base_resolution} {params.zoom_resolutions} '
        '{output.mcool} {input.bg2_files} '
        '>& {log}'

#################################
# Per-celltype bulk mcool (sorted-bulk validation ground truth)
#
# The combined bulk mcool feeds Machima2 training. These per-celltype mcools
# are held-out targets for evaluate.smk's validate_against_sorted_bulk:
# NPC-only and CN-only Bonev tracks act as approximate ground-truth Hi-C of
# each cell-type group, against which Machima's deconvolved per-component
# maps (aggregated by cluster -> celltype) are compared.
#################################
rule bg2_to_mcool_brain_celltype:
    """
    Aggregate the BedGraph2D files of a single celltype group (NPC or CN)
    into a per-celltype multi-resolution mcool. Reuses the already-extracted
    per-GSM bg2.gz files from extract_bonev_track_brain.
    """
    input:
        bg2_files=lambda wc: expand('data/brain/processed/bonev_bg2/{gsm}.bg2.gz',
                                     gsm=config['brain']['query']['celltype_samples'][wc.celltype]),
        chrom_sizes=CHROM_SIZES
    output:
        mcool='data/brain/processed/brain_hic_{celltype}.mcool'
    wildcard_constraints:
        celltype=r'(npc|cn)'
    params:
        base_resolution=min(config['resolutions']),
        zoom_resolutions=','.join(str(r) for r in sorted(set(config['resolutions']),
                                                          reverse=True))
    log:
        'logs/preprocess/bg2_to_mcool_brain_{celltype}.log'
    benchmark:
        'benchmarks/preprocess/bg2_to_mcool_brain_{celltype}.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'bash src/preprocess_bg2_to_mcool.sh '
        '{input.chrom_sizes} {params.base_resolution} {params.zoom_resolutions} '
        '{output.mcool} {input.bg2_files} '
        '>& {log}'

#################################
# Bin ATAC peaks onto fixed genomic bins per chromosome
#################################
rule bin_atac_by_chr:
    """
    Aggregate the ATAC peak counts in multiome_labeled.h5mu onto a uniform
    bin grid (one bin every {resolution} bp), split by chromosome, and
    write one MatrixMarket per chr. The same bin grid is reused for the
    Hi-C side, so Identity T at this resolution holds (n_k == l_k per chr).

    Outputs a directory of {chr}.mtx + chroms.txt + bins.tsv.gz +
    cells.tsv + labels.tsv (when cell_type exists in obs).
    """
    input:
        h5mu='data/{tissue}/processed/multiome_labeled.h5mu'
    output:
        out_dir=directory('data/{tissue}/processed/atac_bins_{resolution}/')
    wildcard_constraints:
        tissue=r'(brain|mtec)',
        resolution=r'\d+'
    params:
        genome=config['genome']
    log:
        'logs/preprocess/bin_atac_by_chr_{tissue}_{resolution}.log'
    benchmark:
        'benchmarks/preprocess/bin_atac_by_chr_{tissue}_{resolution}.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/preprocess_bin_atac_by_chr.py '
        '--input {input.h5mu} --output-dir {output.out_dir} '
        '--genome {params.genome} --resolution {wildcards.resolution} >& {log}'

#################################
# Slice mcool into per-chromosome symmetric bin x bin matrices
#################################
rule bin_hic_by_chr:
    """
    Slice the bulk Hi-C mcool at {resolution} into per-chromosome symmetric
    matrices. Output schema mirrors bin_atac_by_chr (one .mtx per chr +
    chroms.txt). Bin grids are byte-identical to ATAC at the same
    resolution (both use snapatac2.genome.<genome>.chrom_sizes), so
    Identity T can use diag(n_k) directly.
    """
    input:
        mcool='data/{tissue}/processed/{tissue}_hic.mcool'
    output:
        out_dir=directory('data/{tissue}/processed/hic_bins_{resolution}/')
    wildcard_constraints:
        tissue=r'(brain|mtec)',
        resolution=r'\d+'
    log:
        'logs/preprocess/bin_hic_by_chr_{tissue}_{resolution}.log'
    benchmark:
        'benchmarks/preprocess/bin_hic_by_chr_{tissue}_{resolution}.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/preprocess_bin_hic_by_chr.py '
        '--input {input.mcool} --resolution {wildcards.resolution} '
        '--output-dir {output.out_dir} >& {log}'

rule bin_hic_celltype_by_chr:
    """
    Slice per-celltype mcool into per-chromosome symmetric matrices. Same
    bin grid as bin_hic_by_chr; output goes to a celltype-tagged directory
    so it can be loaded independently of the training input.
    """
    input:
        mcool='data/brain/processed/brain_hic_{celltype}.mcool'
    output:
        out_dir=directory('data/brain/processed/hic_bins_{celltype}_{resolution}/')
    wildcard_constraints:
        celltype=r'(npc|cn)',
        resolution=r'\d+'
    log:
        'logs/preprocess/bin_hic_celltype_by_chr_{celltype}_{resolution}.log'
    benchmark:
        'benchmarks/preprocess/bin_hic_celltype_by_chr_{celltype}_{resolution}.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/preprocess_bin_hic_by_chr.py '
        '--input {input.mcool} --resolution {wildcards.resolution} '
        '--output-dir {output.out_dir} >& {log}'

#################################
# Cluster -> celltype annotation via marker genes
#################################
rule annotate_cluster_celltype_brain:
    """
    Score each leiden cluster against NPC and CN marker gene signatures
    (config.celltype_markers.brain). Each cluster is assigned the higher-
    scoring celltype, or 'other' when the score gap is below score_min_diff.
    Outputs a TSV: cluster, celltype, score_npc, score_cn, score_diff.
    """
    input:
        h5mu='data/brain/processed/multiome_labeled.h5mu'
    output:
        tsv='data/brain/processed/cluster_celltype_map.tsv'
    params:
        npc_markers=','.join(config['celltype_markers']['brain']['npc']),
        cn_markers=','.join(config['celltype_markers']['brain']['cn']),
        min_diff=config['score_min_diff']
    log:
        'logs/preprocess/annotate_cluster_celltype_brain.log'
    benchmark:
        'benchmarks/preprocess/annotate_cluster_celltype_brain.txt'
    conda:
        CONDA_PY
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'python src/annotate_cluster_celltype.py '
        '--input {input.h5mu} --output {output.tsv} '
        '--npc-markers {params.npc_markers} --cn-markers {params.cn_markers} '
        '--min-diff {params.min_diff} >& {log}'
