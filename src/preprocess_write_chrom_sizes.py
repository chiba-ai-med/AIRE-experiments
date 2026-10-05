#!/usr/bin/env python
"""Dump <genome>.chrom.sizes from snapatac2.genome.<genome>.

Used to keep the same chromosome name + length list across:
  - misha gdb setup (preprocess_extract_bonev_track.R)
  - cooler load / zoomify  (preprocess_bg2_to_mcool.sh)
  - per-chr binning        (preprocess_bin_atac_by_chr.py / bin_hic_by_chr.py)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import snapatac2 as snap


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--genome", default="mm10")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    g = getattr(snap.genome, args.genome)
    with open(args.output, "w") as f:
        for chrom, length in g.chrom_sizes.items():
            f.write(f"{chrom}\t{length}\n")
    print(f"[write] {args.output}: {len(g.chrom_sizes)} chromosomes")


if __name__ == "__main__":
    main()
