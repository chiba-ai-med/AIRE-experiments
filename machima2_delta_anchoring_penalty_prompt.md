# Machima v1.7.1 — Fix `delta` penalty: anchor to RNA-derived init (`lambda_delta_anchor`)

## Audience

This prompt is for another Claude Code session that will modify the upstream R package **`kokitsuyuzaki/Machima`** (https://github.com/kokitsuyuzaki/Machima). The user is the package maintainer. v1.7.0 (commit `043522e`) implemented `use_shared_background = TRUE` per the spec in `notes/machima2_shared_background_differential_prompt.md`. This prompt fixes a design flaw in that v1.7.0 spec.

## Problem with v1.7.0

The v1.7.0 spec specified the penalty as `lambda_delta · ‖δ_c‖²_F` — a Frobenius-norm penalty that pulls `δ_c` toward **zero**. Implementation in `R/updateDelta.R` faithfully follows this:

```r
denom <- t(T) %*% B %*% (h_c * gd) + L1 + L2 * delta[, c] + 2 * lambda_delta * delta[, c]
```

Empirical result on HiRES (mouse embryo paired scHi-C + scRNA, 2589 cells across 5 cell types) using v1.7.0:

| `lambda_delta` | `‖δ‖_F` (chr1, total of 5 columns) | per-column norm | Outcome |
|---|---|---|---|
| 1 | **4.9e-7** | 0, 0, 0, 0, 0 | δ collapses to zero → R_c = h_c · g_0 · g_0^T (identical per cell type) |
| 0 | 555 | 38, 459, 167, 33, 260 | δ active but **drifts away from cell-type-specific RNA differential**; Hi-C-only loss makes δ_c absorb arbitrary contact signal, not cell-type signature |

In both cases, per-cell-type prediction is **not cell-type-specific**: HiRES specificity (self-Pearson minus mean of others) is in [-0.05, +0.02] for all 5 cell types — same as v1.6.1 baseline.

## Root cause

The penalty `‖δ_c‖²_F` is **unilateral**: it pulls `δ_c` toward 0, while `w_0` (the shared background) is unpenalized. This breaks the intended balance: the optimization always prefers to absorb shared structure into `w_0` (free) rather than `δ_c` (penalized), even when the data signal is genuinely cell-type-specific. With `lambda_delta = 0`, the penalty vanishes but `δ_c` has no inductive bias keeping it tied to cell type `c` — it just becomes one of `J` Hi-C-only basis vectors with no labelling.

The intended design — that `δ_c` reflects the cell-type-c-specific deviation as informed by RNA — is not actualizable with this penalty. The RNA differential is used only to initialise `δ_c`; without an explicit anchoring term, MU updates drift `δ_c` away from the RNA init.

## Fix: anchoring penalty

Replace (or add alongside) the current penalty with an **anchoring** penalty toward the initial `δ_c_init`:

```
new_penalty = lambda_delta_anchor · Σ_c Σ_k ‖δ_c_k − δ_c_init_k‖²_F
```

where `δ_c_init_k` is the value of `delta[[k]][, c]` at iteration 0 (i.e., the RNA-derived init from `pmax(W_RNA[[k]][, c] - rowMeans(W_RNA[[k]]), 1e-5)` or user-supplied `init_delta`).

The gradient of this w.r.t. `δ_c_k` is `2 · lambda_delta_anchor · (δ_c_k - δ_c_init_k)`. Splitting into positive/negative parts for MU:

- Positive part (`+2 · lambda_delta_anchor · δ_c_k`) → denominator
- Negative part (`−2 · lambda_delta_anchor · δ_c_init_k`) → numerator (note minus sign in gradient → plus sign in numerator since we want to pull toward init)

### New MU update for δ_c

```
numer_c = T^T · A · (h_c · (g_0 + δ_c))  +  2 · lambda_delta_anchor · δ_c_init
denom_c = T^T · B · (h_c · (g_0 + δ_c))  +  2 · lambda_delta_anchor · δ_c
          + 2 · lambda_delta · δ_c                 (existing sparsity term, optional)
          + L1 + L2 · δ_c

δ_c ← δ_c · (numer_c / denom_c)^ρ
```

### Behaviour

- `lambda_delta_anchor → ∞`: δ_c is frozen at `δ_c_init` (the RNA differential). Cell-type-specific by construction.
- `lambda_delta_anchor = 0`: no anchor; same as v1.7.0 with `lambda_delta = 0` (δ drifts freely).
- Intermediate values: trade-off between RNA-derived init and Hi-C data fit.

### Default

Recommend `lambda_delta_anchor = 1` as default (moderate anchoring) and keep `lambda_delta = 0` as default (don't shrink δ toward 0 — that was always counter-productive).

## API spec

Add **one** new argument to `Machima2()`:

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `lambda_delta_anchor` | numeric scalar, `≥ 0` or `Inf` | `1` | Strength of anchoring penalty pulling `δ_c` toward its init value `δ_c_init`. `0` = no anchor (v1.7.0 behaviour). `Inf` = freeze `δ_c = δ_c_init`. |

Change default of existing `lambda_delta` from `1` to **`0`** (since the sparsity interpretation was a design error, not user expectation).

`δ_c_init` must be cached at iteration 0 (after `init_delta` resolution) and accessible to `updateDelta`. Suggested implementation: add `delta_init` argument to `.updateDelta` and `.updateDelta_List`, populated from `.initShared` result.

## Files to touch (upstream Machima)

- `R/Machima2.R`: add `lambda_delta_anchor` arg; pass to update loop. Change default `lambda_delta = 0`. Cache `delta_init` from `.initShared` output.
- `R/initShared.R`: also return `delta_init` (currently `list$delta` is both used as init and updated).
- `R/updateDelta.R`: extend signature to accept `delta_init`, `lambda_delta_anchor`. Apply new numerator/denominator terms.
- `man/Machima2.Rd`: document new arg.
- `NEWS.md`: v1.7.1 entry summarizing this fix.
- `tests/testthat/test-shared-background.R`: add test verifying:
  - `lambda_delta_anchor = Inf` keeps `δ` at init values exactly
  - `lambda_delta_anchor = 0` reproduces v1.7.0 behaviour
- `DESCRIPTION`: bump 1.7.0 → 1.7.1

## Validation requirement

Re-run on HiRES embryo 5 cell types after v1.7.1 lands. AIRE-experiments has the data ready and the eval scripts in `/tmp/eval_hires_v170_cross.R`. Required outcome:

- At `lambda_delta_anchor = 1` (default), expect **per-column norms of δ to remain close to init values** (within 2-3× of `‖δ_c_init‖`), not collapse to 0 and not blow up to 500+.
- Cross-Pearson specificity should improve over v1.7.0 baseline ([-0.05, +0.02]) on at least 3/5 cell types into [+0.05, +0.20] range. If fewer than 3 cell types improve, the structural hypothesis (RNA-anchored δ_c → cell-type-specific Hi-C) is falsified and we should abandon this line.

## Pre-merge backup convention

Before modifying `updateDelta.R` and `initShared.R`, copy them to `R/updateDelta.R.preregularized` and `R/initShared.R.preregularized`. Do not delete. Matches existing convention.
