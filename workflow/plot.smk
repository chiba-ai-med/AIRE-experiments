# Plot phase: aggregate metrics across all stage x T_variant runs into a
# single multi-page PDF (background pages explaining the model + result
# pages: Pareto plot, per-chrom heatmap, convergence, sorted-bulk validation).
#
# Driven by src/summarize_machima2.R; the PDF re-renders any time any
# metrics.csv / sorted_validation_brain_*.csv changes.

#################################
# Phase target rule.
#################################
rule plot:
    """Final summary PDF."""
    input:
        f'output/summary_brain_{DEFAULT_RESOLUTION}.pdf',

#################################
# Brain summary PDF
#################################
rule summarize_machima2_brain:
    """
    Render output/summary_brain_<resolution>.pdf with model background +
    cross-stage metrics (Pareto, per-chrom heatmap, convergence, sorted-bulk
    validation).
    """
    input:
        # Self-fit eval metrics (ARI / rel_frob)
        eval_metrics=expand_stage_t(
            'output/eval_brain_{stage}_T{T_variant}_{resolution}/metrics.csv',
            resolution=DEFAULT_RESOLUTION),
        # Held-out per-celltype validation
        sorted_csv=expand_stage_t(
            'output/sorted_validation_brain_{stage}_T{T_variant}_{resolution}.csv',
            resolution=DEFAULT_RESOLUTION),
    output:
        pdf=f'output/summary_brain_{DEFAULT_RESOLUTION}.pdf'
    params:
        eval_glob=f'output/eval_brain_*_{DEFAULT_RESOLUTION}',
        sorted_glob=f'output/sorted_validation_brain_*_{DEFAULT_RESOLUTION}.csv'
    log:
        f'logs/plot/summarize_machima2_brain_{DEFAULT_RESOLUTION}.log'
    benchmark:
        f'benchmarks/plot/summarize_machima2_brain_{DEFAULT_RESOLUTION}.txt'
    conda:
        CONDA_R
    container:
        f"docker://{DOCKER_IMAGE}"
    shell:
        'Rscript src/summarize_machima2.R '
        '"{params.eval_glob}" {output.pdf} "{params.sorted_glob}" >& {log}'
