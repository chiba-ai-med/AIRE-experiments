# Prompt: Investigate H_Sym "background absorbing" component (Part 2/3 follow-up)

## Audience

This prompt is for another Claude Code session that will modify the
**Machima R package** (https://github.com/kokitsuyuzaki/Machima),
v1.0.0 / SHA `c1d1804`. The user is the package maintainer.

## Context

This is the follow-up to `notes/machima2_hsym_dimnames_na_prompt.md`.
That prompt had three parts:

- **Part 1 (done in `c1d1804`)**: never produce NA in H_Sym/H_RNA
  dimnames — fall back to `paste0("comp_", j)`. `assignCelltypeNames.R`
  has been patched.
- **Part 2 (this prompt)**: investigate **why** one component ends up
  as a "background absorbing" component with H_Sym diagonal an order
  of magnitude larger than every other diagonal entry.
- **Part 3 (this prompt)**: document the dimnames-assignment algorithm
  in the function docstring so the meaning of `dimnames(res$H_Sym)` is
  unambiguous.

## Symptom (recap from the original prompt)

In AIRE-experiments / brain 100kb run
(`output/machima2_brain_joint_Tidentity_100000.rds`):

```
H_Sym row/col labels: cluster_4, cluster_1, cluster_0, cluster_6, cluster_3, NA, cluster_2
diagonal values     :   0.184  ,   160.8 ,   139.4 ,    54.4 ,    11.2 , 947.6,    28.7
```

Slot 6 (the previously-NA, now `comp_6` after Part 1 fix) carried
`H_Sym[6,6] = 947` while every other diagonal entry was below 161. Its
off-diagonal entries are uniformly ~200–410 across all pairs. The
component absorbs the bulk variance of the Hi-C signal but is not
assignable to any leiden cluster (max-bipartite-match leaves it
unmatched).

## Hypotheses (in order of likelihood)

1. **Multiplicative-update fixed point**: NMF-style multiplicative
   updates have a known failure mode where one component captures the
   common signal (mean structure) while others specialize. The
   "background" component lives in the null space of the cluster
   structure but has the largest L2 magnitude.

2. **`J` exceeds informative rank**: `J = 7` matches the leiden cluster
   count but the actual rank of the inter-cluster Hi-C contrast might
   be < 7, leaving one component to soak up unstructured variance.
   Test: rerun with `J = 6` and check whether the symptom disappears.

3. **No orthogonality / sparsity pressure on W or H_RNA**: With purely
   non-negative updates and no orthogonality penalty, components are
   free to overlap. Adding a mild `‖W^T W − I‖_F²` orthogonality
   penalty (option `orthW_RNA = TRUE`, already exists in API but its
   strength setting may not be effective) is a standard remedy.

4. **Initialisation seed bias**: `RandomEpi` / `RandomRNA` init draws
   from a distribution that gives one component a head start. Test by
   re-running with `nmf_init_n_restart >= 5` and the
   `init_*` API providing a balanced (e.g., orthogonalized via SVD)
   warm start.

## What to investigate

Run the same brain 100kb input through Machima2 with each of the
following diagnostic configurations, and report which (if any)
eliminate the background-component pattern (defined as: no single
H_Sym diagonal entry is more than 5x the median diagonal entry):

| config | flags                                                              |
|--------|--------------------------------------------------------------------|
| A      | baseline (current default)                                         |
| B      | `J = 6` (one less component than leiden cluster count)             |
| C      | `nmf_init_n_restart = 10`, all else default                        |
| D      | `orthW_RNA = TRUE` (existing API; check that it's wired and effective) |
| E      | `init_W_RNA = svd(stack of X_RNA[k])$u[, 1:J]`-derived warm start  |
| F      | `T_regularization = "frobenius_unit"` (or "low_rank" once available) — does this on its own fix it via H_Sym scale recovery? |

The brain 100kb input is small enough (≈ 2086 cells, J=7, K=22 chroms)
that running all 6 takes a few hours total on the current setup.

## What to fix (if Part 2 yields a clear culprit)

If the diagnostic identifies one of (1)–(4) as the dominant cause:

- **(1) or (2)** → expose a `J` auto-suggestion based on
  `rank-by-cumulative-variance` of `X_RNA + sum_k X_Epi[k]`, with a
  warning when `J` exceeds the suggestion by more than 2.
- **(3)** → verify `orthW_RNA = TRUE` actually adds the
  `‖W^T W − I‖_F²` penalty to the W update (read `R/updateW_RNA2.R`).
  If the flag is silently a no-op, fix it. If it works but is too weak
  by default, consider adding a `lambda_ortho_W` knob.
- **(4)** → bump `nmf_init_n_restart` default from 1 to 3 or 5
  (the cost is bounded by `nmf_init_num_iter * n_restart` on Stage A
  NMF only, so 3-5x is acceptable).

If no single intervention fixes it, document the residual pattern in
the function docstring as a known limitation.

## Part 3: Documentation

Add a `@details` section to `Machima2()`'s Roxygen block describing the
dimnames-assignment algorithm:

```
@details
The dimnames of \code{res$H_Sym} (and rownames of \code{res$H_RNA}) are
assigned by .estimateCelltypes(): for each component j, find the cluster
c that maximizes the correlation of H_RNA[j, ] with the leiden one-hot
indicator of c (across all cells). If the resulting cluster->component
mapping is one-to-one, dimnames[j] is set to that cluster's label.
Otherwise, max_bipartite_match() resolves the mapping; components left
unmatched are labeled "comp_j" (e.g. "comp_6") rather than NA.

Components labeled "comp_*" are typically background-absorbing components:
they capture variance not attributable to any single cluster. Their H_Sym
diagonal entry is often the largest in the matrix. Downstream analyses
that interpret res$H_Sym as a cluster-by-cluster contact propensity should
either drop comp_* rows/cols or treat them as a baseline.
```

## Acceptance criteria

1. A new file `inst/extdata/h_sym_background_diagnostic.md` (or a section
   in vignettes) reports which of A-F eliminates the symptom on the
   brain 100kb fixture, with concrete H_Sym diagonal values for each.
2. `Machima2()` Roxygen `@details` documents the dimnames algorithm
   including the `comp_*` fallback semantics.
3. If a fix is implemented (orthW_RNA wiring, default n_restart bump,
   etc.), unit tests cover it.
4. AIRE-experiments / brain 100kb joint+Tidentity, after pulling the
   fix: `max(diag(H_Sym)) / median(diag(H_Sym)) < 5` for at least one
   of the diagnostic configurations.
5. `R CMD check` clean.

## Out of scope

- Removing the background component from the model (it is informative
  about overall Hi-C scale and removing it can hurt reconstruction).
- Changing the cluster-assignment algorithm itself
  (`max_bipartite_match`); only document and fall-back.
- Refactoring the `init` enum.

## Useful pointers

- Original prompt (Part 1 done, Parts 2-3 carried over here):
  `notes/machima2_hsym_dimnames_na_prompt.md`
- Diagnostic that surfaced the issue:
  `src/diagnose_hsym_and_bonev_diff.R` in AIRE-experiments
- H_Sym summary CSV: `output/diagnostic_brain_100000_hsym_summary.csv`
- H_Sym detail print: `output/diagnostic_brain_100000_hsym_detail.txt`
- Current dimnames logic:
  `R/assignCelltypeNames.R` (`.estimateCelltypes` and the new
  `comp_*` fallback added in `c1d1804`)
- Related upstream work:
  - `notes/machima2_frobenius_unit_listmode_fix_prompt.md` (config F above
    can only be tested cleanly after the frobenius_unit list-mode bug
    is fixed)
  - `notes/machima2_t_regularization_low_rank_prompt.md` (config F's
    `low_rank` option requires this)
