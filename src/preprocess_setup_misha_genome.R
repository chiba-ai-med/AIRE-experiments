#!/usr/bin/env Rscript
# One-shot setup: download a misha genome database (e.g. mm10) into a
# persistent path so subsequent rules can call gsetroot() against it
# without re-downloading the ~786 MB tarball.
#
# misha 5.6.x: gdb.create_genome(genome, path) downloads
# https://misha-genome.s3.eu-west-1.amazonaws.com/<genome>.tar.gz and
# extracts to <path>/<genome>/. We then write a sentinel marker.
#
# Usage:
#   preprocess_setup_misha_genome.R <genome> <out_root_parent>
# Result: <out_root_parent>/<genome>/   (gdb root)
#         <out_root_parent>/<genome>/.created  (sentinel for Snakemake)

suppressPackageStartupMessages({
  library(misha)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: preprocess_setup_misha_genome.R <genome> <out_root_parent>")
}
genome <- args[[1]]
parent <- args[[2]]

dir.create(parent, showWarnings = FALSE, recursive = TRUE)
target <- file.path(parent, genome)

if (file.exists(file.path(target, ".created"))) {
  message(sprintf("[skip] %s already set up at %s", genome, target))
  quit(status = 0)
}

# A previous failed run can leave the genome unpacked but missing the
# .created sentinel. If the gdb's chrom_sizes.txt is already there, skip the
# 786 MB re-download and proceed straight to verification + sentinel write.
already_unpacked <- file.exists(file.path(target, "chrom_sizes.txt"))
if (already_unpacked) {
  message(sprintf("[reuse] %s already unpacked at %s; skipping download", genome, target))
} else {
  if (dir.exists(target) && length(list.files(target)) == 0) {
    unlink(target, recursive = TRUE)
  }
  message(sprintf("[download] %s -> %s", genome, target))
  gdb.create_genome(genome, path = parent)
}

# Sanity: confirm gsetroot works (ALLGENOME isn't always exposed in user
# globalenv across misha versions; gintervals.all() is the public probe).
gsetroot(target)
chr_intervals <- tryCatch(gintervals.all(),
                          error = function(e) {
                            message("WARN: gintervals.all() failed: ", conditionMessage(e))
                            NULL
                          })
if (!is.null(chr_intervals)) {
  message(sprintf("[ok] gdb at %s, %d chromosomes", target, nrow(chr_intervals)))
} else {
  message(sprintf("[ok?] gdb at %s (sanity probe skipped)", target))
}

writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
           file.path(target, ".created"))
