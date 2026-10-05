# Shared constants + module setup. Included from every phase smk
# (preprocess, machima, evaluate, plot) and from the top-level Snakefile.
#
# Loads the project config; do NOT add a configfile: directive in the phase
# smks themselves to avoid double-loading. The Snakefile is the single
# entry point for `configfile:`.

from snakemake.utils import min_version
import glob

min_version("8.10.0")

# Docker image (for --use-singularity mode). Built later from ./Dockerfile.
DOCKER_IMAGE = "koki/aire_workflow:dev"

# Conda environments (for --use-conda mode)
CONDA_R  = "envs/r-machima.yaml"
CONDA_PY = "envs/py-hic.yaml"

# Shared chrom_sizes (genome from config, single source of truth used by
# cooler load/zoomify and per-chr binning).
CHROM_SIZES = f"data/{config['genome']}.chrom.sizes"

# misha genome database root (one-shot ~786 MB download via gdb.create_genome).
# Used by extract_bonev_track to deserialise Bonev StatQuadTreeCached files.
MISHA_GDB = f"data/misha_genomes/{config['genome']}"

# Default (currently the only) resolution. Used as the suffix in per-chr
# bin directories and Machima output paths.
DEFAULT_RESOLUTION = max(config['resolutions'])

# Valid (stage, T_variant) combinations. The matrix is mostly rectangular,
# but jointLowRank requires a learned T (T_regularization=low_rank rejects
# fixT=TRUE in upstream Machima >= 1.1.0), so it only emits Tdense.
def stage_t_pairs():
    pairs = []
    for stage in config['machima2']['stages']:
        if stage == 'jointLowRank':
            t_variants = ['dense']
        else:
            t_variants = config['machima2']['T_variants']
        for tv in t_variants:
            pairs.append((stage, tv))
    return pairs

def expand_stage_t(template, **kw):
    """Expand `template` over valid (stage, T_variant) pairs, filling extras from kw."""
    return [template.format(stage=s, T_variant=tv, **kw) for s, tv in stage_t_pairs()]
