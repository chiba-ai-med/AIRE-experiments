# Machima v1.6.0 — Soft-coupled W between ATAC and Hi-C (`lambda_coupling`)

You are working on the upstream R package **`kokitsuyuzaki/Machima`**. As of v1.5.0, `Machima2()` shares `W_RNA` between the scRNA-side reconstruction and the Hi-C-side reconstruction:

- `X_RNA[[k]] ≈ W_RNA[[k]] · H_RNA`
- `X_Epi[[k]] ≈ (T[[k]] · W_RNA[[k]]) · H_Sym · (T[[k]] · W_RNA[[k]])ᵀ` ( + W_hic block if `J_hic_only > 0`)

Empirical evaluation in AIRE-experiments (brain 100 kb, mTEC 100 kb) showed that this **hard sharing** of `W_RNA` is too restrictive: cell-type columns satisfy ATAC structure but produce near-identical predicted Hi-C contact maps for every cell type (pairwise PC1 correlation ≥ 0.95 across 10 mTEC cell types). The hypothesis the field assumes — that the same low-rank basis describes both ATAC peak structure and Hi-C cell-type-specific contacts — appears not to hold for these datasets at the bin scale tested.

A direct check confirmed the asymmetry: symmetric NMF on bulk mTEC Hi-C alone (no ATAC tie) decomposes into 10 distinct rank-1 patterns with pairwise PC1 correlations spanning −0.70 to +0.88; but loses cell-type labelling. The current `Machima2()` is on the other extreme — labelled but degenerate at the Hi-C side.

This prompt introduces **soft coupling**: allow `W_RNA` and the Hi-C-side basis to differ via a learnable deviation `U`, with a quadratic penalty `lambda_coupling · ‖U‖²_F` controlling how far they may drift. The hard-sharing regime is recovered exactly at `lambda_coupling = Inf` (default), giving full backward compatibility.

## Spec

Three new arguments to `Machima2()`:

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `lambda_coupling` | numeric scalar, `≥ 0` or `Inf` | `Inf` | Coupling strength. `Inf` = hard share (v1.5.0 behaviour). `0` = `W` and `W+U` fully independent (`U` unpenalized). |
| `init_U` | `NULL` or list (length `length(chroms)`) | `NULL` | Optional initial `U[[k]]` per chromosome; each `n_k × J`, non-negative. `NULL` → `U[[k]] = 0`. |
| `fixU` | logical scalar | `is.infinite(lambda_coupling)` | If `TRUE`, `U` is not updated. Auto-defaults to `TRUE` when `lambda_coupling = Inf` (matches v1.5.0 exactly). |

### Extended model

When `lambda_coupling < Inf`, the Hi-C-side basis becomes `W_E[[k]] = W_RNA[[k]] + U[[k]]`:

```
X_RNA[[k]]    ≈ W_RNA[[k]] · H_RNA                            # ATAC unchanged
X_Epi[[k]]    ≈ G_full[[k]] · H_Sym_full · G_full[[k]]ᵀ
                where  G_full[[k]] = [ T[[k]] · (W_RNA[[k]] + U[[k]]) | W_hic[[k]] ]
                and    H_Sym_full  = block_diag(H_Sym, diag(h_hic))

Joint penalty: + lambda_coupling · Σ_k ‖U[[k]]‖²_F
```

`U[[k]]` has the same shape as `W_RNA[[k]]` (`n_k × J`) and is **non-negative**. The W_hic / h_hic mechanism from v1.5.0 is orthogonal and stacks naturally with `U` (both can be active simultaneously).

Result list gains, when `lambda_coupling < Inf`:

- `res$U`: list of length `length(chroms)`, named by chrom; each `n_k × J`, non-negative. Represents the Hi-C-specific deviation of the cell-type basis from the ATAC-derived `W_RNA`.

`res$W_RNA`, `res$H_RNA`, `res$H_Sym`, `res$T`, `res$W_hic`, `res$h_hic` keep their existing shapes and semantics.

### Boundary behaviour

- `lambda_coupling = Inf` and `fixU = TRUE` (defaults): result is **bit-identical to v1.5.0** for any input. `res$U` is **not added** to the result list.
- `lambda_coupling < Inf` and `fixU = TRUE`: `U` is held at `init_U` (or zero), used in Hi-C reconstruction but not updated. Useful as a "fixed prior offset" mode.
- `lambda_coupling < Inf` and `fixU = FALSE`: `U` learnable.
- `lambda_coupling = 0` (or extremely small): `U` is essentially unpenalized; the Hi-C side becomes the bilinear factorisation of `X_Epi` with shape inheriting cell-type column indexing from `W_RNA` only via initialization.

## Math (β-divergence MU; matches existing convention)

Throughout, `R_full[[k]] = (T[[k]]·(W_RNA[[k]] + U[[k]]))·H_Sym·(T[[k]]·(W_RNA[[k]] + U[[k]]))ᵀ + W_hic[[k]]·diag(h_hic)·W_hic[[k]]ᵀ` is the full Hi-C reconstruction at iteration `t`.

Because `W_E = W_RNA + U`, the chain rule gives `∂L_Epi/∂W_RNA = ∂L_Epi/∂U = ∂L_Epi/∂W_E`. The Hi-C numerator/denominator terms are computed once with `W_E` substituted for `W_RNA` and shared between the `W_RNA` and `U` updates.

### `W_RNA[[k]]` update

Combines the ATAC contribution (unchanged) and the Hi-C contribution evaluated at `W_E = W_RNA + U`:

```
num_W[k]   = Pi_RNA[k] · ( X_RNA[k] ⊙ R_RNA[k]^{β-2} ) · H_RNA^T
           + Pi_Epi[k] · ( T[k]^T · (X_Epi[k] ⊙ R_full[k]^{β-2}) · T[k] · (W_RNA[k] + U[k]) · H_Sym )

denom_W[k] = Pi_RNA[k] · R_RNA[k]^{β-1} · H_RNA^T
           + Pi_Epi[k] · ( T[k]^T · R_full[k]^{β-1} · T[k] · (W_RNA[k] + U[k]) · H_Sym )
           + L1_W_RNA + L2_W_RNA · W_RNA[k]

W_RNA[k] ← W_RNA[k] · ( num_W[k] / denom_W[k] )^{ρ(β)}
```

(`R_RNA[k] = W_RNA[k] · H_RNA` is the ATAC reconstruction; β-divergence MU with relaxation `ρ(β)` matches the existing `R/updateW_RNA2.R` convention.)

Note that `W_RNA` continues to feel Hi-C loss gradient through the `(W_RNA + U)` substitution. This is intentional and matches the chain rule. The L2 penalty on `U` (below) is what biases the model to put deviation into `U` only when ATAC alone cannot explain Hi-C.

### `U[[k]]` update (new)

Hi-C loss only, plus quadratic penalty `lambda_coupling · ‖U[[k]]‖²_F`. The penalty contributes `2 · lambda_coupling · U[k]` to the denominator (gradient `d‖U‖²/dU = 2U`, positive part since `U ≥ 0`).

```
num_U[k]   = Pi_Epi[k] · ( T[k]^T · (X_Epi[k] ⊙ R_full[k]^{β-2}) · T[k] · (W_RNA[k] + U[k]) · H_Sym )

denom_U[k] = Pi_Epi[k] · ( T[k]^T · R_full[k]^{β-1} · T[k] · (W_RNA[k] + U[k]) · H_Sym )
           + L1_U + L2_U · U[k]
           + 2 · lambda_coupling · U[k]

U[k] ← U[k] · ( num_U[k] / denom_U[k] )^{ρ(β)}
```

`L1_U`, `L2_U` are **reused from `L1_W_RNA`, `L2_W_RNA`** — do not introduce new parameters. (The coupling penalty is the dominant regularizer on `U`.)

### Other updates

`H_RNA`, `H_Sym`, `T`, `W_hic`, `h_hic` are updated by their existing rules, with one substitution: wherever the Hi-C reconstruction or Hi-C gradient currently uses `W_RNA`, use `W_RNA + U` instead. Specifically:

- `R/updateH_Sym.R` / `R/updateH_Sym_diag.R`: `G_shared = T · (W_RNA + U)` instead of `T · W_RNA`. `R_full` definition unchanged form-wise.
- `R/updateT.R`: gradient uses `W_RNA + U` as the basis being mapped.
- `R/updateW_hic.R`, `R/updateH_hic.R`: `R_full` now includes the `U` contribution to the shared block. No change to the `W_hic` / `h_hic` MU form.

### `lambda_coupling = Inf` short-circuit

When `lambda_coupling = Inf`:
- `fixU` defaults to `TRUE` (set in `.checkMachima2`)
- `U` initialized to zero (matrix of zeros) and never updated
- `W_RNA + U = W_RNA` throughout, so all updates reduce to v1.5.0 exactly
- `res$U` is **omitted** from the result list

Implementation: a single `if (is.infinite(lambda_coupling)) skip U-related code paths` guard at the top of the iteration loop is sufficient. Behaviour must be bit-identical to v1.5.0 with the same seed.

## Tasks

1. **`R/Machima2.R`** — argument list: add `lambda_coupling = Inf`, `init_U = NULL`, `fixU = NULL` (sentinel for auto-default) near `J_hic_only`. Pass through `.checkMachima2`, `.initMachima2`, and the iteration loop. After every existing reconstruction site that previously computed `T %*% W_RNA` for Hi-C, switch to a helper `.G_shared(W_RNA, U, T)` returning `T %*% (W_RNA + U)`. Attach `U` to the result only when `lambda_coupling < Inf`.

2. **`R/checkMachima2.R`** — validate:
   - `is.numeric(lambda_coupling) && length == 1 && (is.infinite(lambda_coupling) || lambda_coupling >= 0)`
   - if `!is.null(init_U)`: list of length `length(chroms)`, each `n_k × J`, non-negative
   - if `is.null(fixU)`: set `fixU <- is.infinite(lambda_coupling)` (auto-default)
   - `is.logical(fixU) && length == 1`
   - error if `fixU = TRUE && lambda_coupling < Inf && is.null(init_U)` is allowed (U = 0 is a valid fixed state)

3. **`R/initMachima2.R`** — when `lambda_coupling < Inf`:
   - if `init_U` supplied, use it
   - else: per chrom, `U[[k]] <- matrix(0, n_k, J)` (zero init: starts at hard-share regime, deviates only as data warrants)
   - reuse the same RNG path so `set.seed()` reproducibility holds
   - when `lambda_coupling = Inf`: `U` is `NULL` (signal to update loop to skip U code paths)

4. **`R/updateU.R`** (new file) — implement the MU rule above. Signature mirrors `.updateW_RNA2` (list mode, per-chrom Reduce). Reuse computed `R_full[[k]]` from the helper.

5. **`R/updateW_RNA2.R`** — modify the Hi-C contribution to use `W_RNA[k] + U[k]` instead of `W_RNA[k]` when `U` is not NULL. No change when `U = NULL`. The ATAC contribution is unchanged.

6. **`R/updateH_Sym.R`** and **`R/updateH_Sym_diag.R`** — same substitution `W_RNA → W_RNA + U` in Hi-C gradient computation.

7. **`R/updateT.R`** — same substitution.

8. **`R/updateW_hic.R`** and **`R/updateH_hic.R`** — `R_full` already includes the W_hic block; just ensure the shared block uses `W_RNA + U`.

9. **`R/reconstructEpi.R`** (or wherever `.reconstructEpi_single` lives, added in v1.5.0) — update the `G_shared` calculation to `T_k %*% (W_RNA_k + U_k)` when `U_k` is non-NULL, else `T_k %*% W_RNA_k`.

10. **`man/Machima2.Rd`** — add `@param lambda_coupling`, `@param init_U`, `@param fixU`. Document new return-value entry `U` under `@return`.

11. **Tests** — add `tests/testthat/test-Machima2-soft-coupled-W.R`:
    - `lambda_coupling = Inf` (default) reproduces v1.5.0 exactly (bit-identical with same seed).
    - `lambda_coupling < Inf, fixU = TRUE, init_U = NULL` reproduces v1.5.0 exactly (U fixed at 0).
    - `lambda_coupling < Inf, fixU = FALSE` produces non-NULL `res$U` with shape matching `W_RNA` per chrom.
    - As `lambda_coupling` decreases, `RecError` on `X_Epi` decreases monotonically (more flexibility → better fit).
    - `init_U` shape validation: wrong shape → error.
    - `lambda_coupling = 0` produces a valid non-negative `U` (no NaN / Inf).

12. **NEWS** — `inst/NEWS`:
    ```
    VERSION 1.6.0
    ------------------------
       o Added lambda_coupling parameter to Machima2() for soft coupling
         between the ATAC-side W_RNA and the Hi-C-side basis. When finite,
         introduces a learnable deviation U with W_E := W_RNA + U for the
         Hi-C reconstruction, regularized by lambda_coupling * ||U||_F^2.
         Default Inf reproduces v1.5.0 hard-sharing behaviour exactly.
       o Added init_U and fixU arguments mirroring the W_hic_init / fixW_hic
         pattern from v1.5.0.
    ```

13. **`DESCRIPTION`** — bump `Version: 1.6.0`.

## Out of scope

- Do **not** add `lambda_coupling` support to the asymmetric `Machima()` function. This change is `Machima2()`-only.
- Do **not** introduce per-cell-type `lambda_coupling[c]` in this version. Single scalar is sufficient for the initial deployment; per-column variants can be added later if empirically motivated.
- Do **not** add separate L1/L2 parameters for `U` (`L1_U`, `L2_U`); reuse `L1_W_RNA` and `L2_W_RNA`.
- Do **not** apply T to U separately — `U` lives in the same `n_k` space as `W_RNA[[k]]`, and gets the same `T_k` transformation via `T_k %*% (W_RNA + U)`.
- Do **not** make `U` symmetric across `W_RNA` columns (e.g. swap U's column 1 and 2 — the column ordering of `U` must match `W_RNA`).

## Acceptance

- `R CMD check` clean.
- All existing tests pass with defaults (bit-exact identical to v1.5.0).
- New `test-Machima2-soft-coupled-W.R` passes.
- On a small synthetic example (e.g. 50 bin × 30 cell ATAC, 50 × 50 Hi-C, J=4, J_hic_only=0), running with `lambda_coupling = 1.0` yields `tail(RecError, 1) < tail(RecError_v1.5.0, 1)` strictly (more parameters → better fit), and `‖U‖_F > 0`.

Commit subject: short (Soft-coupled W via U deviation); body: explain motivation (hard W sharing forces cell-type columns to compromise between ATAC and Hi-C objectives; empirically observed degeneracy in AIRE-experiments brain 100 kb and mTEC 100 kb where predicted per-cell-type Hi-C maps have PC1 correlation ≥ 0.95 across all cell types; soft coupling via lambda_coupling-regularized deviation `U` allows the Hi-C basis to drift from `W_RNA` proportionally to data evidence while preserving cell-type column ordering through the regularization-anchored initialization at U=0).
