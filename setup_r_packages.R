#!/usr/bin/env Rscript

#################################
# Install R packages for AIRE-experiments
#################################
# Run AFTER conda env (workflow/envs/r-machima.yaml) is created.
# Conda already provides Matrix, ggplot2, fs (with libuv), remotes, devtools,
# rcppeigen, rann, rspectra (compiled deps for nnTensor's GitHub variant).
# Here we add only:
#   - nnTensor : GitHub (rikenbit/nnTensor); the GitHub branch ships
#                graph-regularized NMF (L_graph_U/V, lambda_graph_*,
#                algorithm "TV-GNMF") not present in CRAN 1.3.0
#   - Machima  : GitHub-only (kokitsuyuzaki/Machima)
#   - misha    : Tanay genome-DB; needed to deserialise Bonev 2017 Hi-C tracks
#
# Vicus is pulled in transitively as a rikenbit/nnTensor dependency
# (graph.method="Vicus" requires it).
#
# Seurat / Signac / SummarizedExperiment / SingleCellExperiment / MuData /
# rhdf5 are intentionally NOT installed -- no R rule in this project uses
# them (Python handles all multiome I/O).

cat("Installing R packages for AIRE-experiments...\n\n")

cran_packages <- character()  # all R packages now via conda or GitHub

# GitHub specs with optional minimum version pin. When a min_version is given
# and the installed version is older, the package is reinstalled from HEAD of
# the default branch (which is expected to satisfy the pin). Bump min_version
# here whenever a new upstream release adds a feature this repo depends on.
github_packages <- list(
  list(spec = "rikenbit/nnTensor",     name = "nnTensor", min_version = NULL),
  list(spec = "kokitsuyuzaki/Machima", name = "Machima",  min_version = "1.7.1"),
  list(spec = "tanaylab/misha",        name = "misha",    min_version = NULL)
)

cat("=== CRAN ===\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org/", dependencies = TRUE)
  } else {
    cat(pkg, "already installed.\n")
  }
}

cat("\n=== GitHub ===\n")
for (pkg in github_packages) {
  installed   <- requireNamespace(pkg$name, quietly = TRUE)
  needs_inst  <- !installed
  if (installed && !is.null(pkg$min_version)) {
    have <- packageVersion(pkg$name)
    if (have < pkg$min_version) {
      cat(pkg$name, "installed", as.character(have),
          "<", pkg$min_version, "(required) -- reinstalling\n")
      needs_inst <- TRUE
    }
  }
  if (needs_inst) {
    cat("Installing", pkg$spec, "...\n")
    remotes::install_github(pkg$spec, dependencies = TRUE, upgrade = "never")
  } else {
    cat(pkg$name, "already installed (",
        as.character(packageVersion(pkg$name)), ").\n")
  }
}

cat("\n=== Verify ===\n")
all_pkgs <- c(cran_packages, vapply(github_packages, function(p) p$name, character(1)))
ok <- TRUE
for (pkg in all_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("OK  ", pkg, as.character(packageVersion(pkg)), "\n")
  } else {
    cat("FAIL", pkg, "\n")
    ok <- FALSE
  }
}

if (!ok) {
  quit(status = 1)
}
cat("\nAll packages installed.\n")
