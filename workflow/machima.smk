# Machima phase: Stage A (NMF) + Stage B (Machima2) per stage x T_variant.
#
# Depends on per-chromosome ATAC + Hi-C bin matrices produced by the
# preprocess phase. One rule (run_machima2) parameterised by:
#   stage     : joint | transferFlog | transferKraw | supervisedHfix
#               | supervisedGNMF | supervisedWinit
#   T_variant : identity | dense
#
# Stage semantics live in src/run_machima2.R; see the comment block at the
# top of that script for what each stage does at the algorithm level.
#
# Wildcards / config / shared constants come from common.smk (Snakefile-loaded).

#################################
# Phase target rule.
#################################
# Stage-x-T_variant matrix isn't fully rectangular -- see stage_t_pairs() in
# common.smk. jointLowRank only emits Tdense (fixT=TRUE conflicts with
# T_regularization=low_rank in upstream Machima >= 1.1.0).
rule machima:
    """All trained Machima2 .rds outputs (brain, J=7, every valid stage x T)."""
    input:
        expand_stage_t(
            'output/machima2_brain_{stage}_T{T_variant}_{resolution}.rds',
            resolution=DEFAULT_RESOLUTION),

#################################
# Machima2 deconvolution
#################################
rule run_machima2:
    """
    Call run_machima2.R on per-chromosome ATAC + Hi-C lists.

    Wildcards:
      {tissue}    = brain | mtec
      {stage}     = joint | transferFlog | transferKraw
                  | supervisedHfix | supervisedGNMF | supervisedWinit
                  | jointLowRank
      {T_variant} = identity | dense (jointLowRank: dense only)
      {resolution} in bp

    Saves a single .rds with the Machima2 result list (W_RNA, H_RNA, H_Sym,
    T, RecError, RelChange, .meta).
    """
    input:
        atac_dir='data/{tissue}/processed/atac_bins_{resolution}/',
        hic_dir='data/{tissue}/processed/hic_bins_{resolution}/'
    output:
        rds='output/machima2_{tissue}_{stage}_T{T_variant}_{resolution}.rds'
    wildcard_constraints:
        tissue=r'(brain|mtec)',
        stage=r'(joint|transferFlog|transferKraw|supervisedHfix|supervisedGNMF|supervisedWinit|jointLowRank)',
        T_variant=r'(identity|dense)',
        resolution=r'\d+'
    params:
        J=config['machima2']['J'],
        num_iter=config['machima2']['num_iter'],
        transfer_num_iter=config['machima2']['transfer_num_iter'],
        nmf_num_iter=config['machima2']['nmf_num_iter'],
        nmf_n_restart=config['machima2']['nmf_n_restart'],
        gnmf_lambda_V=config['machima2']['gnmf_lambda_V'],
        T_regularization=config['machima2']['T_regularization'],
        lambda_T=config['machima2']['lambda_T'],
        T_rank=config['machima2']['T_rank'],
        H_Sym_structure=config['machima2']['H_Sym_structure'],
        lambda_balance=config['machima2']['lambda_balance']
    log:
        'logs/machima/run_machima2_{tissue}_{stage}_T{T_variant}_{resolution}.log'
    benchmark:
        'benchmarks/machima/run_machima2_{tissue}_{stage}_T{T_variant}_{resolution}.txt'
    conda:
        CONDA_R
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'Rscript src/run_machima2.R '
        '{input.atac_dir} {input.hic_dir} {wildcards.stage} {wildcards.T_variant} '
        '{output.rds} {params.J} {params.num_iter} {params.transfer_num_iter} '
        '{params.nmf_num_iter} {params.nmf_n_restart} {params.gnmf_lambda_V} '
        '{params.T_regularization} {params.lambda_T} {params.T_rank} '
        '{params.H_Sym_structure} {params.lambda_balance} '
        '>& {log}'
