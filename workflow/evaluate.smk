# Evaluate phase: per-Machima-run metrics + sorted-bulk validation.
#
#   evaluate_machima2          : ARI vs leiden + per-chrom Frobenius
#                                reconstruction (training-fit only)
#   validate_against_sorted_bulk: per-celltype reconstruction compared
#                                against held-out NPC-only / CN-only Bonev
#                                bulk Hi-C (Pearson, cosine, scale-invariant
#                                rel_frob)

#################################
# Phase target rule.
#################################
rule evaluate:
    """All metrics needed before plot.smk can render."""
    input:
        # Self-fit eval (ARI + per-chrom Frobenius)
        expand_stage_t(
            'output/eval_brain_{stage}_T{T_variant}_{resolution}/metrics.csv',
            resolution=DEFAULT_RESOLUTION),
        # Held-out NPC/CN bulk validation
        expand_stage_t(
            'output/sorted_validation_brain_{stage}_T{T_variant}_{resolution}.csv',
            resolution=DEFAULT_RESOLUTION),

#################################
# Self-fit evaluation
#################################
rule evaluate_machima2:
    """
    Convergence trace + cell-type alignment (H_RNA argmax vs leiden ARI)
    + per-chromosome reconstruction Frobenius error. Writes metrics.csv,
    convergence_*.pdf, confusion_matrix.csv, reconstruction_per_chr.csv/pdf
    under output/eval_<tissue>_<stage>_T<variant>_<resolution>/.
    """
    input:
        rds='output/machima2_{tissue}_{stage}_T{T_variant}_{resolution}.rds',
        atac_dir='data/{tissue}/processed/atac_bins_{resolution}/',
        hic_dir='data/{tissue}/processed/hic_bins_{resolution}/'
    output:
        metrics='output/eval_{tissue}_{stage}_T{T_variant}_{resolution}/metrics.csv'
    wildcard_constraints:
        tissue=r'(brain|mtec)',
        stage=r'(joint|transferFlog|transferKraw|supervisedHfix|supervisedGNMF|supervisedWinit|jointLowRank)',
        T_variant=r'(identity|dense)',
        resolution=r'\d+'
    params:
        out_dir='output/eval_{tissue}_{stage}_T{T_variant}_{resolution}'
    log:
        'logs/evaluate/evaluate_machima2_{tissue}_{stage}_T{T_variant}_{resolution}.log'
    benchmark:
        'benchmarks/evaluate/evaluate_machima2_{tissue}_{stage}_T{T_variant}_{resolution}.txt'
    conda:
        CONDA_R
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'Rscript src/evaluate_machima2.R '
        '{input.rds} {input.atac_dir} {input.hic_dir} {params.out_dir} >& {log}'

#################################
# Held-out per-celltype validation (sorted-bulk ground truth)
#################################
rule validate_against_sorted_bulk:
    """
    Compare Machima2's per-celltype reconstruction to held-out per-celltype
    Bonev bulk Hi-C. For each component j, infer its dominant leiden cluster
    from H_RNA, map cluster -> {npc, cn, other} via cluster_celltype_map.tsv,
    aggregate the components in each celltype group via H_Sym masking:

        X_hat_g[k] = (T[k] W[k]) M_g H_Sym M_g^T (T[k] W[k])^T

    Compares per-chrom: Pearson correlation, cosine similarity, and a
    scale-invariant relative-Frobenius (each side normalised by its own
    Frobenius norm before subtraction).
    """
    input:
        rds='output/machima2_brain_{stage}_T{T_variant}_{resolution}.rds',
        atac_dir='data/brain/processed/atac_bins_{resolution}/',
        npc_dir='data/brain/processed/hic_bins_npc_{resolution}/',
        cn_dir='data/brain/processed/hic_bins_cn_{resolution}/',
        cluster_map='data/brain/processed/cluster_celltype_map.tsv'
    output:
        csv='output/sorted_validation_brain_{stage}_T{T_variant}_{resolution}.csv'
    wildcard_constraints:
        stage=r'(joint|transferFlog|transferKraw|supervisedHfix|supervisedGNMF|supervisedWinit|jointLowRank)',
        T_variant=r'(identity|dense)',
        resolution=r'\d+'
    log:
        'logs/evaluate/validate_against_sorted_bulk_brain_{stage}_T{T_variant}_{resolution}.log'
    benchmark:
        'benchmarks/evaluate/validate_against_sorted_bulk_brain_{stage}_T{T_variant}_{resolution}.txt'
    conda:
        CONDA_R
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'Rscript src/validate_against_sorted_bulk.R '
        '{input.rds} {input.atac_dir} {input.npc_dir} {input.cn_dir} '
        '{input.cluster_map} {output.csv} >& {log}'
