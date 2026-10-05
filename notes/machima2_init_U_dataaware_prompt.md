# Machima v1.6.1 — Data-aware default initialization for `U` (soft-coupling)

You are working on the upstream R package **`kokitsuyuzaki/Machima`**. v1.6.0 introduced `lambda_coupling` and the deviation matrix `U` for soft-coupled W between ATAC/RNA and Hi-C bases:

- `X_RNA[[k]] ≈ W_RNA[[k]] · H_RNA`
- `X_Epi[[k]] ≈ (T · (W_RNA + U)[[k]]) · H_Sym · (T · (W_RNA + U)[[k]])ᵀ`
- Penalty: `lambda_coupling · ‖U‖²_F`

Empirical evaluation in AIRE-experiments mTEC 100 kb revealed a critical issue: **when `init_U = NULL` (default) and `lambda_coupling < Inf`, `U` is initialized to zero and never escapes**. The multiplicative update for `U` has the form

```
U ← U · (num / denom)^{ρ(β)}
```

so `U_init = 0 ⟹ U^{(t)} = 0` for all `t`. This silently degrades the soft-coupling feature to `lambda_coupling = Inf` regardless of the user-supplied value. Verified directly: a λ scan across `{1000, 100, 10, 1}` gave bit-identical results (all with `‖U‖_F = 0`) on the mTEC dataset under the default init.

The fix: when `lambda_coupling < Inf` and `init_U = NULL`, initialize `U` to small **non-zero** values matched to the typical magnitude of `W_RNA`, so the MU has a starting point from which it can grow or shrink as data dictates.

## Spec

The behaviour of the existing arguments is preserved (no API change at the user level). Only the internal default initialization in `R/initMachima2.R` changes.

### Current (v1.6.0) default

```r
# in initMachima2 List branch, when J_hic_only > 0 ... wait no, that was W_hic
# For U:
if (lambda_coupling < Inf && is.null(init_U)) {
    U <- lapply(X_RNA, function(x) matrix(0, nrow(x), J))
}
```

(Zero initialization. Causes MU dead-zone.)

### New (v1.6.1) default

```r
if (lambda_coupling < Inf && is.null(init_U)) {
    U <- lapply(seq_along(X_RNA), function(k) {
        # Scale: 1/100 of W_RNA's max-init magnitude.
        # W_RNA is conventionally initialized via runif(.,0,1) at start of
        # NMF inner loop; we use a small fraction of that scale for U so
        # the initial Hi-C reconstruction contribution from (W_RNA + U)
        # is W_RNA-dominated (preserving the v1.5.0-like starting point)
        # while U has non-zero values that the MU can update.
        n_k <- nrow(X_RNA[[k]])
        # Match the W_RNA init scale used by .initMachima2_List/.initMachima2_Matrix
        # (uniform [0,1] divided by 100 for U):
        matrix(runif(n_k * J, 1e-5, 1e-2), n_k, J)
    })
}
```

(Uniform random `[1e-5, 1e-2]`. Non-zero, much smaller than typical `W_RNA` values, so the starting reconstruction is essentially the v1.5.0 W_RNA-only contribution plus a tiny perturbation.)

The corresponding `.initMachima2_Matrix` branch (single-matrix mode) gets the analogous one-shot init.

### Why these specific values

The mTEC empirical workaround used `init_U ~ runif(1e-5, 1e-4)` and yielded successful MU progression (`||U||_F / ||W||_F` ratios from 0.0002 at λ=1000 to 200+ at λ=0.01). Values in `[1e-5, 1e-2]` cover this range; we anchor the upper bound at `1e-2` to be safe for datasets where `W_RNA` itself converges to slightly larger magnitudes (e.g., RNA-binned log-normalised data has values up to ~5, `W_RNA` may converge to magnitudes ~0.1).

Using `runif(., 1e-5, 1e-2)` rather than a constant ensures the starting point breaks symmetry across the J columns of `U` (otherwise all columns of `U` would be updated identically in the first iteration and the model would be effectively rank-1 from the start).

### Reproducibility

The init uses `runif`, so it consumes the RNG state. Place this draw **after** the `W_RNA`/`H_RNA`/`H_Sym` initialization draws so that:
- `set.seed(s)` at the call site reproduces existing behaviour for `lambda_coupling = Inf`
- For `lambda_coupling < Inf`, the same seed produces the same init for both old (zero) and new (non-zero) defaults *plus* the user-supplied `init_U` path is unchanged.

## Tasks

1. **`R/initMachima2.R`** — in both `.initMachima2_Matrix` and `.initMachima2_List` branches, locate the block that allocates `U` when `lambda_coupling < Inf && is.null(init_U)`. Replace `matrix(0, ...)` (List) or `matrix(0, ...)` (Matrix) with `matrix(runif(n_k * J, 1e-5, 1e-2), n_k, J)`. The RNG call must come after the existing `W_RNA`/`H_RNA`/`H_Sym` initialization so that the seed semantics for the `lambda_coupling = Inf` path are unaffected.

2. **`R/Machima2.R`** — update the `@param init_U` line of the roxygen doc to reflect that the auto-default is non-zero uniform random.

3. **`man/Machima2.Rd`** — re-rendered via `devtools::document()` from the updated roxygen.

4. **Tests** — add `tests/testthat/test-Machima2-U-init.R`:
   - `lambda_coupling = Inf` (default) reproduces v1.5.0 exactly (bit-identical with same seed).
   - `lambda_coupling = 1`, `init_U = NULL`: `res$U[[k]]` after iteration 1 contains entries of order ~1e-5 to ~1e-2 (the post-MU values, which depend on the data; just assert > 0).
   - `lambda_coupling = 1`, two runs with same seed produce identical `res$U`.
   - `lambda_coupling = 1`, two runs with different seeds produce *different* `res$U` (verifies RNG is consumed correctly).
   - `lambda_coupling = 1`, supplying `init_U` explicitly overrides the auto-default.

5. **NEWS** — `inst/NEWS`:
    ```
    VERSION 1.6.1
    ------------------------
       o Fixed default initialization of U when lambda_coupling < Inf.
         The previous default (U_init = 0) interacted with the
         multiplicative update to keep U at zero throughout training,
         silently disabling the soft-coupling feature. The new default
         draws U_init from runif(., 1e-5, 1e-2), small enough to leave
         the starting Hi-C reconstruction essentially unchanged from
         v1.5.0 but non-zero so the MU can move U as data warrants.
       o init_U argument behaviour unchanged when supplied explicitly.
       o lambda_coupling = Inf path bit-identical to v1.5.0.
    ```

6. **`DESCRIPTION`** — bump `Version: 1.6.1`.

## Out of scope

- Do **not** change the math of the `U` update rule. The fix is purely an initialization scale change.
- Do **not** change `lambda_coupling`, `fixU`, or any other v1.6.0 arguments.
- Do **not** add data-dependent scaling (e.g., scale to `mean(X_Epi)^0.5 / J`). The fixed range `[1e-5, 1e-2]` is sufficient for the empirical regime we tested; data-dependent scaling can be added in a future release if needed.

## Acceptance

- `R CMD check` clean.
- All existing tests pass.
- New `test-Machima2-U-init.R` passes.
- On the mTEC J=10 setup, `lambda_coupling ∈ {1000, 100, 10, 1}` with `init_U = NULL` produces visibly different `‖U‖_F` values across the four runs (not all zero).

Commit subject: short (Fix U init dead-zone); body: explain that v1.6.0's `U_init = 0` interacted with the multiplicative update to silently disable soft coupling, and that the fix is a small non-zero uniform initialization.
