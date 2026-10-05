# Prompt: Fix NA in H_Sym dimnames + investigate "background absorbing" component in Machima2

## Audience

This prompt is for another Claude Code session that will modify the
**Machima2 R package** (https://github.com/kokitsuyuzaki/Machima).
The user is the package maintainer.

## Context

`Machima2()` writes per-component labels into the dimnames of `H_Sym`
(and the rownames of `H_RNA`). These labels are intended to map each
component j (1..J) to the dominant cluster of its top-loaded cells, so
downstream code can read `dimnames(res$H_Sym)[[1]][j]` to learn which
cluster each component "represents".

## Symptom

In AIRE-experiments / brain 100 kb run (`output/machima2_brain_joint_Tidentity_100000.rds`):

```r
> dimnames(res$H_Sym)
[[1]]
[1] "cluster_4" "cluster_1" "cluster_0" "cluster_6" "cluster_3"  NA
[7] "cluster_2"

[[2]]
[1] "cluster_4" "cluster_1" "cluster_0" "cluster_6" "cluster_3"  NA
[7] "cluster_2"
```

**One of the J=7 components has dimname `NA`.** The user provided
non-NA cluster labels for all 2086 cells (verified: `cells.tsv` and
`labels.tsv` 100% match, no NA in the cluster column). So the NA must
be generated inside Machima2's own labeling logic.

## Why this matters

Looking at the actual H_Sym values for that component (slot 6):

```
H_Sym row/col labels: cluster_4, cluster_1, cluster_0, cluster_6, cluster_3, NA, cluster_2

diagonal values     :       0.184 ,    160.8,    139.4,     54.4 ,    11.2 , 947.6,    28.7
```

The NA-labeled component has H_Sym diagonal **947** while every other
diagonal is below 161. Its row/column entries are uniformly large
(~200-410). This component dominates H_Sym across all pairs.

Interpretation: it looks like a **"background absorbing component"**
that captured the bulk variance of the Hi-C signal, and Machima2's own
component-to-cluster majority vote could not assign it to any specific
cluster (so it labeled it NA).

This is not strictly a crash, but:

- Downstream `dimnames`-based analysis sees NA and may silently drop it
- The component biases per-celltype reconstructions in unexpected ways
- It is not clear whether this is a *bug* or a *modeling choice*

## What to investigate / fix

### Part 1 (definite fix): never produce NA dimnames

Whatever logic populates H_Sym / H_RNA dimnames should fall back to
`paste0("comp_", j)` when no dominant cluster can be assigned.
Acceptance: `any(is.na(dimnames(res$H_Sym)[[1]]))` returns FALSE for
every output of `Machima2()`.

### Part 2 (investigate): why does one component end up "background"

Possible causes, in order of likelihood:

1. **Random init bias**: with `joint` stage on a moderately balanced
   input, the optimizer may consolidate variance into a single
   high-magnitude component while letting others specialize. This
   is a well-known NMF failure mode.
2. **Insufficient sparsity / orthogonality** in the model loss — adding
   a small `||W^T W - I||_F^2` orthogonality penalty often fixes it
3. **Component count J too high** for the actual number of celltypes
   in the data — the extra component absorbs noise

Determine which applies and document. If (1) or (2), consider adding
an option (off by default) to inject mild orthogonality regularization
on W or H_RNA to discourage background components.

### Part 3 (investigate): cluster assignment for the NA component

How does Machima2 internally assign the dimname? Is it majority vote
of which.max(H_RNA[j, ]) cells, or something else? Document the
algorithm in the function docstring so the meaning of dimnames is
unambiguous.

## Out of scope

- Changing the labeling algorithm itself (just make it never return NA)
- Changing model output values (only dimnames)
- Removing the "background component" (if it is unavoidable, just label
  it cleanly)

## Acceptance test

The user will re-run AIRE-experiments / brain 100 kb (joint stage,
both T variants) and check:

```r
res <- readRDS("output/machima2_brain_joint_Tidentity_100000.rds")
stopifnot(!any(is.na(dimnames(res$H_Sym)[[1]])))
stopifnot(!any(is.na(rownames(res$H_RNA))))
```

If Part 2 leads to an orthogonality option: with that option enabled,
no single H_Sym diagonal entry should be more than ~5x the median
H_Sym diagonal magnitude (rough heuristic for "no single component
absorbs all variance").

## Useful pointers

- Diagnostic that surfaced the issue: `src/diagnose_hsym_and_bonev_diff.R`
  in AIRE-experiments
- H_Sym summary CSV: `output/diagnostic_brain_100000_hsym_summary.csv`
- H_Sym detail print (with the NA visible):
  `output/diagnostic_brain_100000_hsym_detail.txt`
- Related issue (separate prompt):
  `notes/machima2_dense_t_regularization_prompt.md`
