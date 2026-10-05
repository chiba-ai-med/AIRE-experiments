#!/usr/bin/env Rscript
# Extract a Bonev/Tanay misha 2D track tarball into a BedGraph2D (bg2.gz)
# at the project's finest bin resolution.
#
# Bonev et al. 2017 (GSE96107) ships per-GSM Hi-C as a directory of
# StatQuadTreeCached binary files (one per chromosome pair) inside a
# `<sample>.track/` directory, packaged as `<GSM>_<sample>.tar.gz`. These
# are misha 2D tracks; only the misha R package (tanaylab/misha) can
# deserialise them.
#
# Pipeline:
#   1. Untar into a fresh misha gdb (with mm10 chrom_sizes pre-written)
#   2. Move the extracted .track dir under <gdb>/tracks/
#   3. gdb.reload() to register
#   4. Per chromosome pair, gextract at iterator=c(res, res) -> data.frame
#      (chrom1, start1, end1, chrom2, start2, end2, value)
#   5. Append rows to <output>.bg2.gz (BedGraph2D, gzipped)
#
# The bg2.gz is consumed by `cooler load -f bg2 mm10:res in.bg2 out.cool`.
# Iterator size = finest bin resolution = 10 kb, since coarser bins are
# obtained for free via `cooler zoomify` later.
#
# Usage:
#   preprocess_extract_bonev_track.R <tarball> <chrom_sizes> <resolution> <out_bg2_gz>

suppressPackageStartupMessages({
  library(misha)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: preprocess_extract_bonev_track.R <tarball> <misha_gdb_root> <resolution> <out_bg2_gz>")
}
tarball       <- normalizePath(args[[1]])
misha_root    <- normalizePath(args[[2]])  # e.g. data/misha_genomes/mm10
resolution    <- as.integer(args[[3]])
# misha::gsetroot() changes the working directory, which silently breaks any
# subsequent relative path. Make the output absolute before doing anything
# that may chdir.
out_bg2_gz    <- normalizePath(args[[4]], mustWork = FALSE)
if (!startsWith(out_bg2_gz, "/")) {
  out_bg2_gz <- file.path(getwd(), out_bg2_gz)
}

stopifnot(file.exists(tarball), dir.exists(misha_root), resolution > 0)

# Stop write.table from emitting "1e+05" instead of "100000" for bin coords —
# cooler's bg2 parser tolerates scientific notation in some versions but the
# explicit-integer form is unambiguous.
options(scipen = 999)

dir.create(dirname(out_bg2_gz), showWarnings = FALSE, recursive = TRUE)

#################################
# 1. Activate the pre-built misha gdb (created once by setup_misha_genome rule)
#################################
gsetroot(misha_root)
message(sprintf("[gdb] %s", misha_root))

#################################
# 2. Untar the Bonev tarball; move the .track dir into the gdb
#################################
untar_dir <- tempfile("bonev_untar_")
dir.create(untar_dir, recursive = TRUE)
message(sprintf("[untar] %s -> %s", tarball, untar_dir))
untar(tarball, exdir = untar_dir)

track_dirs <- list.files(untar_dir, pattern = "\\.track$",
                         full.names = TRUE, recursive = TRUE,
                         include.dirs = TRUE)
track_dirs <- track_dirs[file.info(track_dirs)$isdir %in% TRUE]
if (length(track_dirs) == 0) {
  stop(sprintf("no *.track directory found inside %s", tarball))
}
track_src   <- track_dirs[[1]]
track_name  <- sub("\\.track$", "", basename(track_src))
track_dst   <- file.path(misha_root, "tracks", paste0(track_name, ".track"))
# If a stale copy from a previous run exists, drop it.
if (dir.exists(track_dst)) unlink(track_dst, recursive = TRUE)
message(sprintf("[move] %s -> %s", track_src, track_dst))
dir.create(dirname(track_dst), showWarnings = FALSE, recursive = TRUE)
# file.rename() fails across filesystems (e.g. /tmp tmpfs -> data/ disk).
# Use system mv which falls back to copy+unlink automatically.
mv_status <- system2("mv", c(shQuote(track_src), shQuote(track_dst)))
if (mv_status != 0L) {
  stop(sprintf("mv failed: %s -> %s (status=%d)",
               track_src, track_dst, mv_status))
}

gdb.reload()
if (!(track_name %in% gtrack.ls())) {
  stop(sprintf("track %s not visible after gdb.reload(); known tracks: %s",
               track_name, paste(gtrack.ls(), collapse = ", ")))
}
message(sprintf("[track] registered as '%s'", track_name))

#################################
# 3. Per chromosome pair, gextract -> append to bg2.gz
#################################
# misha's <gdb>/chrom_sizes.txt uses bare names ("1", "2", ...) while
# gintervals.all() (and downstream cooler) use the "chr" prefix. Use
# gintervals.all() as the canonical source so the names line up everywhere.
chr_intervals <- gintervals.all()
# misha returns chrom as a factor whose levels include non-canonical contigs
# (chrM, *_random, chrUn_*, ...). Subsetting drops the rows but NOT the unused
# levels, so `chroms[[j]]` later returns a factor whose integer level code can
# exceed length(chrom_size), giving "subscript out of bounds" on chrX/chrY.
# Coerce to plain character before doing any indexing.
chr_intervals <- chr_intervals[grepl("^chr[0-9XY]+$", as.character(chr_intervals$chrom)), ]
chroms <- as.character(chr_intervals$chrom)
chrom_size <- setNames(chr_intervals$end, chroms)
message(sprintf("[chroms] %d canonical (max=%s)",
                length(chroms), names(which.max(chrom_size))))

out_conn <- gzfile(out_bg2_gz, open = "w")

n_total_rows <- 0L
for (i in seq_along(chroms)) {
  c1 <- chroms[[i]]
  for (j in seq(i, length(chroms))) {
    c2 <- chroms[[j]]
    intervals <- gintervals.2d(c1, 0, chrom_size[[c1]],
                               c2, 0, chrom_size[[c2]])
    df <- tryCatch(
      gextract(track_name, intervals = intervals,
               iterator = c(resolution, resolution),
               colnames = "value"),
      error = function(e) {
        message(sprintf("  [skip] %s-%s: %s", c1, c2, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(df) || nrow(df) == 0) next

    # Some misha versions emit chrom1/chrom2/start1/start2/end1/end2/value;
    # others use chrom/start/end1/start2/end2 -- normalise.
    needed <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "value")
    missing <- setdiff(needed, colnames(df))
    if (length(missing) > 0) {
      stop(sprintf("gextract did not return expected columns; missing: %s; got: %s",
                   paste(missing, collapse = ", "),
                   paste(colnames(df), collapse = ", ")))
    }
    df <- df[, needed, drop = FALSE]
    # Bonev/Tanay misha tracks store float-valued normalised contact intensities
    # and emit NA for empty bins. Cooler's bg2 loader rejects NA + casts count
    # to int32 by default; drop empties and keep float format for the value
    # column (cooler load gets --field count:dtype=float on the consumer side).
    df <- df[!is.na(df$value), , drop = FALSE]
    if (nrow(df) == 0) next
    df$start1 <- as.integer(df$start1)
    df$end1   <- as.integer(df$end1)
    df$start2 <- as.integer(df$start2)
    df$end2   <- as.integer(df$end2)
    write.table(df,
                file = out_conn, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    n_total_rows <- n_total_rows + nrow(df)
    if ((i + j) %% 5 == 0) {
      message(sprintf("  %s-%s: %d rows (cumulative %d)",
                      c1, c2, nrow(df), n_total_rows))
    }
  }
}

flush(out_conn)
close(out_conn)
message(sprintf("[done] %s: %d total rows", out_bg2_gz, n_total_rows))
