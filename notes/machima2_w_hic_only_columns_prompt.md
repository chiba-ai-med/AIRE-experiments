# Machima v1.5.0 — Hi-C-only basis columns (`J_hic_only`)

You are working on the upstream R package **`kokitsuyuzaki/Machima`**. `Machima2()` currently shares `W_RNA` between the scRNA reconstruction and the Hi-C reconstruction:

- `X_RNA[[k]]  ≈  W_RNA[[k]] · H_RNA`
- `X_Epi[[k]]  ≈  G[[k]] · H_Sym · G[[k]]ᵀ`,  `G[[k]] = T[[k]] · W_RNA[[k]]`

Empirical evaluation in AIRE-experiments brain 100 kb showed that **even after fixing rank-1 collapse** (via diagonal masking + log1p preprocessing, which spread `H_Sym` diagonal from `[200, 1980]` to `[200, 413]`), compartment-level recovery against sorted-bulk Hi-C remained at PC1 correlation ≈ 0.04 (≈ noise floor) with diff PC1 ≈ -0.01. The bottleneck is **structural**: the cell-type basis `W_RNA` is forced to simultaneously explain ATAC peak structure (per-cell, per-bin) and Hi-C chromatin contact structure (per-bin-pair) with the same column vectors. These two objectives have non-overlapping signal in many genomic regions, so `W_RNA` cannot capture Hi-C-specific structure (e.g. compartment boundaries that don't align with ATAC peaks).

This prompt adds **Hi-C-only basis columns** that participate only in the Hi-C reconstruction, leaving `W_RNA` to specialise on the cell-type basis.

## Spec

New parameter: **`J_hic_only ∈ Z_≥0`** (default `0L`, backward-compatible — `0L` reduces to the v1.4.x model exactly).

### Extended model

When `J_hic_only > 0`:

```
G_full[[k]]   = [ T[[k]] · W_RNA[[k]] | W_hic[[k]] ]        # l_k × (J + J_hic_only)
H_Sym_full    = block_diag( H_Sym , diag(h_hic) )           # (J+J_hic_only) × (J+J_hic_only)
X_Epi[[k]]    ≈ G_full[[k]] · H_Sym_full · G_full[[k]]ᵀ
              = G_shared[[k]] H_Sym G_shared[[k]]ᵀ
              + W_hic[[k]] diag(h_hic) W_hic[[k]]ᵀ
```

where:

- `W_hic[[k]]` is `l_k × J_hic_only`, **non-negative**, learned only via Hi-C loss (no ATAC gradient). One matrix per chromosome (consistent with `W_RNA` list mode).
- `h_hic ∈ ℝ^{J_hic_only}_+` is shared across chromosomes (consistent with `H_Sym`).
- The hic-only block of `H_Sym_full` is **always diagonal** (regardless of `H_Sym_structure`) — `W_hic` columns have no cell-type semantic, so off-diagonal interactions between hic-only components have no biological interpretation.
- Cross-block entries (shared × hic-only) are **forced to zero**. The two reconstructions add cleanly: the shared block carries cell-type structure; the hic-only block absorbs Hi-C-specific residual. This preserves interpretability of `W_RNA` and `H_Sym` as cell-type quantities.
- **`T` is not applied to `W_hic`.** `W_hic` lives directly in Hi-C bin space (`l_k`). Even when `T` is dense (learned), `W_hic` is independent of the ATAC↔Hi-C bin mapping.

### API

`Machima2()` argument list gains:

```r
J_hic_only = 0L,                # number of Hi-C-only basis columns
W_hic_init = NULL,              # optional list, length(chroms); each l_k × J_hic_only non-negative
fixW_hic   = FALSE              # if TRUE, do not update W_hic (analogous to fixW_RNA)
```

Result list gains, when `J_hic_only > 0`:

- `res$W_hic`: list of length `length(chroms)`, named by chrom; each `l_k × J_hic_only`.
- `res$h_hic`: numeric vector length `J_hic_only` (the hic-only diagonal of `H_Sym_full`).

`res$W_RNA`, `res$H_RNA`, `res$H_Sym`, `res$T` keep their existing shapes (J-based). `H_Sym` is **always** the J×J shared block — never the extended `(J+J_hic_only)` matrix. This keeps downstream code that reads `res$H_Sym` working unchanged.

## Math (β-divergence MU; matches existing convention)

Let `R_full[[k]] = G_full[[k]] · H_Sym_full · G_full[[k]]ᵀ` denote the full reconstruction at iteration `t`.

### `W_hic[[k]]` update (new)

For β-divergence with relaxation exponent `ρ(β)`:

```
num_W_hic   = (X_Epi[[k]] ⊙ R_full[[k]]^{β-2}) · W_hic[[k]] · diag(h_hic)
denom_W_hic = R_full[[k]]^{β-1}            · W_hic[[k]] · diag(h_hic) + L1 + L2 · W_hic[[k]]
W_hic[[k]] ← W_hic[[k]] · ( num_W_hic / denom_W_hic )^{ρ(β)}
```

For β = 2 this collapses to:

```
num_W_hic   = X_Epi[[k]] · W_hic[[k]] · diag(h_hic)
denom_W_hic = R_full[[k]] · W_hic[[k]] · diag(h_hic) + L1 + L2 · W_hic[[k]]
```

### `h_hic` update (new) — proper diagonal MU (same form as v1.3.0 `H_Sym` diagonal)

For each `i ∈ {1, …, J_hic_only}`:

```
num_i   = Σ_k  w_hic_i^{(k)ᵀ} ( X_Epi[[k]] ⊙ R_full[[k]]^{β-2} ) w_hic_i^{(k)}
denom_i = Σ_k  w_hic_i^{(k)ᵀ}              R_full[[k]]^{β-1}     w_hic_i^{(k)}
h_hic_i ← h_hic_i · ( num_i / denom_i )^{ρ(β)}
```

L1/L2 on `h_hic` translate as in standard NMF (`denom_i ← denom_i + L1_h_hic`, `denom_i ← denom_i + L2_h_hic · h_hic_i`). Use the **same `L1_H_Sym`, `L2_H_Sym` parameters** as the shared block — do not introduce new parameters for the hic-only block.

### `W_RNA[[k]]`, `H_RNA`, `H_Sym`, `T` updates

Mathematically unchanged in form. The only change is that all reconstructions of `X_Epi[[k]]` must use `R_full[[k]]` (which now includes the hic-only contribution) **in the denominator**. The numerator still uses `X_Epi[[k]]`. The hic-only contribution in `R_full` acts as a fixed background offset for the shared-block updates — exactly the standard MU treatment of additive offsets.

When `J_hic_only == 0`, `W_hic` is empty, `h_hic` is length 0, `R_full` reduces to the current `R_shared`, and all updates are bit-identical to v1.4.x.

## Tasks

1. **`R/Machima2.R`** — argument list: add `J_hic_only = 0L`, `W_hic_init = NULL`, `fixW_hic = FALSE` near `fixW_RNA`. Pass through `.checkMachima2`, `.initMachima2`, and the iteration loop. After every existing reconstruction site that previously computed `G %*% H_Sym %*% t(G)` for `X_Epi`, switch to a helper `.reconstructEpi(W_RNA, T, W_hic, H_Sym, h_hic)` that returns `R_full[[k]]`. In the result list, attach `W_hic` and `h_hic` only when `J_hic_only > 0`.

2. **`R/checkMachima2.R`** — validate:
   - `is.numeric(J_hic_only) && length == 1 && J_hic_only == as.integer(J_hic_only) && J_hic_only >= 0`
   - if `!is.null(W_hic_init)`: list of length `length(chroms)`, each `l_k × J_hic_only`, non-negative
   - `is.logical(fixW_hic) && length == 1`
   - error if `fixW_hic = TRUE && is.null(W_hic_init)` (no initialization to fix to)

3. **`R/initMachima2.R`** — when `J_hic_only > 0`:
   - if `W_hic_init` supplied, use it
   - else: per chrom, `W_hic[[k]] <- matrix(runif(l_k * J_hic_only, 0.1, 1.0), l_k, J_hic_only)` (matching the existing `W_RNA` random-init scale; reuse the same RNG path so `set.seed()` reproducibility holds)
   - initialize `h_hic <- rep(1, J_hic_only)` (matching the existing `H_Sym` diag init convention)

4. **`R/updateW_hic.R`** (new file) — implement the MU rule above. Signature mirrors `.updateW_RNA2` (list mode, per-chrom Reduce). Use `R_full[[k]]` from the helper.

5. **`R/updateH_hic.R`** (new file) — implement diagonal MU on `h_hic`. Reuse the constrained-MU pattern from `R/updateH_Sym.R`'s diagonal branch (added in v1.3.0).

6. **`R/updateW_RNA2.R`** — change reconstruction-of-Epi computation to call `.reconstructEpi`. The `denom1` term that involves `X_Epi` reconstruction now uses `R_full` instead of `G H_Sym Gᵀ`. No other change.

7. **`R/updateH_Sym.R`** — same: replace inline reconstruction with `.reconstructEpi`. The shared-block update (whether `H_Sym_structure="symmetric"` or `"diagonal"`) acts only on the J×J block and uses `G_shared` (= `T %*% W_RNA`), but the **denominator** picks up the hic-only contribution via `R_full`. Operationally:
   - `num_H_Sym[i,j] = Σ_k g_shared_i^{(k)ᵀ} (X_Epi ⊙ R_full^{β-2}) g_shared_j^{(k)}`
   - `denom_H_Sym[i,j] = Σ_k g_shared_i^{(k)ᵀ}              R_full^{β-1}     g_shared_j^{(k)}`

8. **`R/updateT.R`** — same replacement (use `R_full` in denominator). T receives gradient only from the shared block (W_hic does not flow through T).

9. **`R/reconstructEpi.R`** (new helper, internal) —
   ```r
   .reconstructEpi <- function(W_RNA_k, T_k, W_hic_k, H_Sym, h_hic) {
     G_shared <- T_k %*% W_RNA_k
     R <- G_shared %*% H_Sym %*% t(G_shared)
     if (length(h_hic) > 0L) {
       R <- R + W_hic_k %*% (h_hic * t(W_hic_k))   # diag(h_hic) %*% t(W_hic_k) but cheaper
     }
     R
   }
   ```

10. **`man/Machima2.Rd`** — add `@param J_hic_only`, `@param W_hic_init`, `@param fixW_hic`. Document new return-value entries (`W_hic`, `h_hic`) under `@return`.

11. **Tests** — add `tests/testthat/test-Machima2-J-hic-only.R`:
    - `J_hic_only = 0L` (default) reproduces v1.4.x exactly (bit-identical with same seed).
    - `J_hic_only = 2L` produces non-negative `W_hic[[k]]` with `ncol == 2` per chrom, length-2 `h_hic`.
    - `J_hic_only = 2L` reduces `RecError` on `X_Epi` compared to `J_hic_only = 0L` (more parameters → strictly better fit).
    - `H_Sym` returned in result is still `J × J` (not extended).
    - `fixW_hic = TRUE` with supplied `W_hic_init` keeps `W_hic` exactly equal to init across iterations.
    - `W_hic_init` shape validation: wrong `J_hic_only` → error; wrong `l_k` → error; negative entries → error.

12. **NEWS** — `inst/NEWS`:
    ```
    VERSION 1.5.0
    ------------------------
       o Added J_hic_only parameter to Machima2() for Hi-C-only basis
         columns (W_hic). When J_hic_only > 0, the Hi-C reconstruction
         becomes G_full H_full G_full^T where G_full = [T W_RNA | W_hic]
         and H_full is block-diagonal (cross blocks zero). The hic-only
         block is always diagonal; W_hic learns Hi-C-specific structure
         independent of the ATAC cell-type basis.
       o Added W_hic_init and fixW_hic arguments mirroring W_RNA_init
         and fixW_RNA.
       o Default J_hic_only = 0L reproduces v1.4.x behaviour exactly.
    ```

13. **`DESCRIPTION`** — bump `Version: 1.5.0`.

## Out of scope

- Do **not** add `J_hic_only` support to the asymmetric `Machima()` function. This change is `Machima2()`-only.
- Do **not** add cross-block entries to `H_Sym_full`. Cross blocks remain zero by spec.
- Do **not** apply `T` to `W_hic`. Even with `fixT = FALSE`, `W_hic` lives directly in Hi-C bin space.
- Do **not** add label-derived initialisation for `W_hic` (the `label` argument applies only to `W_RNA`, by design — hic-only columns have no cell-type semantic).
- Do **not** introduce `J_hic_only`-specific L1/L2 parameters; reuse `L1_H_Sym`, `L2_H_Sym`, `L1_W_RNA`, `L2_W_RNA` (the latter applied to `W_hic` as well).

## Acceptance

- `R CMD check` clean.
- All existing tests pass with `J_hic_only = 0L` (default) — bit-exact identical to v1.4.x.
- New `test-Machima2-J-hic-only.R` passes.
- On a small synthetic example, `J_hic_only = 2` strictly reduces `tail(RecError, 1)` compared to `J_hic_only = 0`.

Commit subject: short (Hi-C-only basis columns); body: explain motivation (shared `W_RNA` cannot capture Hi-C-specific structure not aligned with ATAC peaks; verified empirically on AIRE-experiments brain 100 kb where compartment recovery stalls at PC1 corr ≈ 0.04 even after rank-1 collapse fix).
