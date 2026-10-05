# Machima v1.7.0 — Shared background + cell-type differential decomposition (`g_0`, `delta_c`)

## Audience

This prompt is for another Claude Code session that will modify the upstream R package **`kokitsuyuzaki/Machima`** (https://github.com/kokitsuyuzaki/Machima). The user is the package maintainer. As of v1.6.1, `Machima2()` has soft-coupled `W` via `lambda_coupling` and data-aware `U` init.

## Background

`Machima2()` was validated on three benchmarks in AIRE-experiments:

| Dataset | Setup | Result |
|---|---|---|
| Brain 100kb | J=3 (npc/cn/other), mixed bulk Hi-C, Bonev sorted truth | Cross-Pearson specificity ±0.04 (failed) |
| Bonev synthetic 3-stage | J=3, in silico-mixed mES+NPC+CN, known 1/3 each | Specificity ±0.05 (failed); H_diag balanced but predicted Hi-C is bulk-mimicking |
| **HiRES J=5** | **5 cell types (blood, ExE endoderm, early neurons, radial glias, neural ectoderm), paired scRNA+scHi-C (real, same cells)** | **Specificity ±0.04 (failed)** |

Across all three, the symptom is identical:

- `H_diag` is balanced (no one-component collapse) at `lambda_coupling = 1`
- BUT each per-cell-type prediction `R_c = h_c · g_c · g_c^T` (with `g_c = T·(W[:,c] + U[:,c])`) is **bulk-mimicking**: the 5 (or 3) predicted Hi-C maps are essentially identical, none of them cell-type-specific against the corresponding sorted/aggregated truth
- `||U||_F / ||W||_F ≈ 5-12x`: U dominates W, so RNA-derived basis has little effect on per-cell-type Hi-C prediction
- H_Sym = symmetric (off-diagonal allowed) gives partial improvement (off-diagonal absorbs ~75% of mass) but predicted slots still collapse to 1-2 differential axes, not J independent cell types

**Differential evaluation** (compute `Δ_c = R_c − mean_c(R_c)` and compare to `Δ_truth_c = truth_c − mean_c(truth_c)`) reveals that on HiRES, Machima2 captures **one** rank-1 differential signal (the dominant cell-type axis, e.g., "early-neurons vs blood") and distributes it across all 5 cell-type slots — but cannot split it into 5 truly independent cell-type signals.

## Root cause

The current model

```
X_Epi[[k]] ≈ G_k · H_Sym · G_kᵀ,   G_k = T_k · (W_k + U_k)
R_c        := h_c · g_c · g_cᵀ  (per-cell-type prediction, where g_c = c-th col of G_k)
```

cannot prevent each `g_c` from absorbing the shared "bulk-like" structure (compartments, distance decay, etc.) of `X_Epi`. With H_Sym = diagonal, the J rank-1 outer products `g_c g_cᵀ` independently add up to reconstruct `X_Epi`, and the optimization preferentially makes each `g_c` look like a slight variation of the bulk — hence bulk-mimicking. Adding U (`lambda_coupling`) gives the model flexibility, but U has no cell-type structure, so it absorbs the shared component into a soup distributed over all J slots.

The structural fix is to **explicitly factor the Hi-C model into a shared background `g_0` and small cell-type-specific deviations `δ_c`**, so that cell-type information is forced to live in the `δ_c` and cannot be drowned by the shared structure.

## Proposed model (v1.7.0)

```
X_Epi[[k]] ≈ h_0 · g_0_k · g_0_kᵀ                       (shared bulk background)
           + Σ_{c=1}^{J} h_c · (g_0_k + δ_c_k) · (g_0_k + δ_c_k)ᵀ   (per-cell-type)

g_0_k = T_k · w_0_k    (length n_k)              — shared basis, no cell-type label
δ_c_k = T_k · ε_c_k    (n_k × 1 each, J columns) — cell-type-specific deviation

penalty: + lambda_delta · Σ_{c,k} ‖δ_c_k‖²
```

`w_0_k` and `ε_c_k` are the new free parameters. The RNA-side reconstruction stays the same shape:

```
X_RNA[[k]] ≈ W_RNA_k · H_RNA       (unchanged; W_RNA columns are cell-type-specific RNA profiles)
```

The link between RNA and Hi-C now flows specifically into `δ_c`:

```
ε_c_k  is initialized from / optionally tied to  W_RNA_k[:, c] − mean_c(W_RNA_k[:, c])
```

i.e., the RNA differential profile per cell type drives the Hi-C differential basis. `w_0_k` is purely a Hi-C nuisance variable (background), with no RNA tie.

### Per-cell-type prediction

```
R_c := h_c · (g_0_k + δ_c_k) · (g_0_k + δ_c_k)ᵀ
     = h_c · ( g_0_k · g_0_kᵀ                      # shared
              + g_0_k · δ_c_kᵀ + δ_c_k · g_0_kᵀ    # cross terms
              + δ_c_k · δ_c_kᵀ )                   # pure differential

Δ_c := R_c − (R̂_avg)  where R̂_avg = (Σ_c h_c R_c) / Σ_c h_c
                       collapses to a function of δ_c minus mean_c(δ_c).
```

The user-facing "predicted cell-type-c Hi-C" is `R_c`. The user-facing "predicted differential Hi-C" is `Δ_c`. **Both must be returned by Machima2** because the latter is the metric on which cell-type specificity is meaningful.

### What this fixes

- Bulk-mimicking: `g_0` is the explicit container for shared structure. The optimizer no longer needs to distribute the bulk across the J cell-type slots — the slots are forced to be deviations.
- J independence: each `δ_c` is independently penalized (its own `‖·‖²` term), so they cannot collapse to a common axis.
- RNA-Hi-C coupling: `δ_c` is tied to RNA differential through initialization (and optionally L2-soft-tied during updates).
- Backward compatibility: at `lambda_delta = Inf` and a degenerate `g_0 = 0` initialisation, the model collapses to current v1.6.1 behaviour (no shared component).

## API spec

Five new arguments to `Machima2()`:

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `use_shared_background` | logical scalar | `FALSE` | Master switch. `FALSE` = v1.6.1 behaviour, ignore other new args. `TRUE` = activate `g_0` + `δ` decomposition. |
| `init_g0` | `NULL` or list (length `length(chroms)`) | `NULL` | Optional `g_0_k` init; each `n_k`-vector, non-negative. `NULL` → derive from first symNMF singular vector of `X_Epi[[k]]`. |
| `init_delta` | `NULL` or list (length `length(chroms)`) | `NULL` | Optional `ε_c_k` init; each `n_k × J`, non-negative. `NULL` → derive from `W_RNA_k[:, c] − mean_c W_RNA_k[:, c]` clipped to non-negative (since model parameters are NN). |
| `lambda_delta` | numeric scalar, `≥ 0` or `Inf` | `1` | Penalty `lambda_delta · ‖δ_c_k‖²`. `Inf` forces `δ = 0` (bulk-only); `0` removes constraint. |
| `fix_g0` | logical scalar | `FALSE` | If `TRUE`, `w_0_k` not updated. |

Result list gains, when `use_shared_background = TRUE`:

- `res$w_0`: list of length `length(chroms)`, named by chrom; each `n_k`-vector, non-negative. The shared background basis.
- `res$delta`: list, named by chrom; each `n_k × J`, non-negative. The cell-type-specific deviations.
- `res$predict_celltype_R(c, chrom)`: convenience function returning `R_c` for one cell type.
- `res$predict_celltype_delta(c, chrom)`: convenience function returning `Δ_c` (the differential map; this is what cell-type-specificity metrics should use).
- The existing `res$W_RNA`, `res$H_RNA`, `res$H_Sym`, `res$T`, `res$U` remain unchanged in shape and semantics. `res$U` is **not** active when `use_shared_background = TRUE` (mutually exclusive with `lambda_coupling`); document this.

## Math (β-divergence MU; matches existing convention)

`H_Sym` should be kept diagonal (`H_Sym_structure = "diagonal"`) in the v1.7.0 path; the off-diagonal role is now played by the cross-terms `g_0 δ_cᵀ + δ_c g_0ᵀ` which the model fits explicitly.

Let

```
R_full[[k]] = h_0 · g_0_k · g_0_kᵀ + Σ_c h_c · (g_0_k + δ_c_k)(g_0_k + δ_c_k)ᵀ
```

at iteration t. Define

```
A_k := X_Epi[k] ⊙ R_full[k]^{β-2}    (data-tied numerator factor)
B_k := R_full[k]^{β-1}                (denominator factor)
```

### `w_0_k` update (new)

The gradient of the Hi-C loss w.r.t. `w_0_k` collects all 1+J contributions:

```
∂R_full/∂w_0 = 2 · h_0 · g_0_k · w_0_kᵀ
            + 2 · Σ_c h_c · (g_0_k + δ_c_k) · w_0_kᵀ
```

Plug into β-divergence MU. Since `g_0_k = T_k · w_0_k`, factor out `T_k`:

```
num_w0[k]   = T_kᵀ · A_k · ( h_0 · g_0_k + Σ_c h_c · (g_0_k + δ_c_k) )
denom_w0[k] = T_kᵀ · B_k · ( h_0 · g_0_k + Σ_c h_c · (g_0_k + δ_c_k) )
            + L1_g0 + L2_g0 · w_0_k

w_0_k ← w_0_k · ( num_w0[k] / denom_w0[k] )^{ρ(β)}
```

(`ρ(β)` is the existing relaxation exponent in `R/updateW_RNA2.R`.)

### `ε_c_k` update (new)

The gradient w.r.t. column `c`:

```
∂R_full/∂ε_c = 2 · h_c · (g_0_k + δ_c_k) · ε_c_kᵀ      (only the c-th outer product contributes)

num_eps[k, :, c]   = T_kᵀ · A_k · h_c · (g_0_k + δ_c_k)
denom_eps[k, :, c] = T_kᵀ · B_k · h_c · (g_0_k + δ_c_k)
                   + L1_delta + L2_delta · ε_c_k
                   + 2 · lambda_delta · ε_c_k

ε_c_k ← ε_c_k · ( num_eps[k, :, c] / denom_eps[k, :, c] )^{ρ(β)}
```

### `H_Sym` update (modified)

`h_0` is treated as a fresh diagonal entry alongside `h_1, ..., h_J`. The update is the same form as the existing v1.6.1 H_Sym update, but now over `J+1` entries with the `(J+1)-th` basis being `g_0_k` directly (no `δ`).

### `W_RNA` update (unchanged for RNA side, modified for Hi-C side)

RNA side: same as v1.6.1 (no Hi-C loss flows into `W_RNA` once `use_shared_background = TRUE`).

This decoupling is the key simplification: with `use_shared_background = TRUE`, Hi-C reconstruction is built from `w_0` and `δ`, neither of which is `W_RNA`. The link is only at initialisation, so `lambda_coupling` is not needed in this mode.

If the user wants `W_RNA` to feel Hi-C loss too, they should use v1.6.1 (`use_shared_background = FALSE` + `lambda_coupling < Inf`). The two modes are exclusive.

## Boundary behaviour

- `use_shared_background = FALSE` (default): exactly v1.6.1.
- `use_shared_background = TRUE, lambda_delta = Inf`: `δ_c = 0`, so all cell types share the same `g_0`; per-cell-type predictions become `R_c = (h_0 + h_c) · g_0 · g_0ᵀ` (cell types differ only by scalar `h_c`). This is the "bulk reconstruction only" limit and useful as a sanity check.
- `use_shared_background = TRUE, lambda_delta = 0`: `δ` unpenalized; should approximately recover bulk reconstruction quality matching v1.6.1 if data permits.

## Validation requirement

Before merging, the v1.7.0 path must be validated on:

1. **HiRES embryo 5 cell types** (data in AIRE-experiments at `data/hires/processed/`). Required: `R_c` cross-Pearson specificity for at least 3/5 cell types in [+0.05, +0.20] AND `Δ_c` specificity matching truth direction for at least 3/5 cell types. The existing v1.6.1 fails on both axes.
2. **Bonev synthetic 3-stage** (`data/brain/processed/rna_bins_bonev3_100000/` + Hi-C). Required: each predicted celltype's self-Pearson must exceed its mean off-Pearson by ≥ 0.05.
3. **Brain ncx (npc/cn)** (existing data). Should not regress (specificity ≥ current ±0.04).

The validation harness in AIRE-experiments is in `/tmp/eval_hires_cross.R` and `/tmp/eval_hires_differential.R`. The maintainer should re-run them after upstream code lands.

## Files to touch (upstream Machima)

- `R/Machima2.R`: add args, dispatch on `use_shared_background`.
- `R/initShared.R` (new): build `w_0_k`, `ε_c_k` initial values.
- `R/updateSharedG.R` (new): MU update for `w_0_k`.
- `R/updateDelta.R` (new): MU update for `ε_c_k`.
- `R/updateH_Sym2.R`: extend to `J+1` entries when `use_shared_background = TRUE`.
- `R/assembleR.R` (helper): compute `R_full` from `w_0`, `δ`, `H_Sym`, `T` for both training (loss eval) and prediction (`predict_celltype_R/delta`).
- `man/Machima2.Rd`: document new args.
- `NEWS.md`: v1.7.0 entry summarizing this prompt.
- `tests/testthat/test-shared-background.R` (new): unit test for backward compat at default, and basic shape test at `use_shared_background = TRUE`.
- `DESCRIPTION`: bump version 1.6.1 → 1.7.0.

## Pre-merge backup convention

Before modifying any update files, copy the v1.6.1 versions with `.preregularized` suffix (e.g., `R/updateW_RNA2.R` → `R/updateW_RNA2.R.preregularized`). Do not delete. This matches the existing convention noted in AIRE-experiments memory.
