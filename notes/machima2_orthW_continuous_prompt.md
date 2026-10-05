# Machima v1.4.0 — Continuous orthogonality strength for `orthW_RNA`

You are working on the upstream R package **`kokitsuyuzaki/Machima`**. `Machima()` and `Machima2()` currently have a binary `orthW_RNA = TRUE/FALSE` flag controlling whether `W_RNA` columns are pushed toward orthogonality. Empirical evaluation in AIRE-experiments showed:

- `orthW_RNA = FALSE`: 1 W column dominates (becomes "background absorber"), other columns drowned out
- `orthW_RNA = TRUE`: balances column magnitudes but reconstruction error increases 4× (over-constrained)

A **continuous strength parameter** is needed to interpolate between these extremes.

## Spec

Replace the binary toggle with a continuous parameter:

- **`lambda_orthW ∈ [0, 1]`** (new; default 0)
  - 0 = standard NMF denominator (current `orthW_RNA = FALSE`)
  - 1 = orthogonal NMF denominator (current `orthW_RNA = TRUE`)
  - intermediate = weighted blend

The existing `orthW_RNA` parameter is **kept for backward compatibility** with deprecation warning:
- `orthW_RNA = TRUE` → automatically sets `lambda_orthW = 1` (with warning)
- `orthW_RNA = FALSE` (default) → `lambda_orthW` retains user-supplied value (or default 0)

## Math

Current MU denominator branch in `R/updateW_RNA2.R` (and `R/updateW_RNA.R` for `Machima()`):

```r
if (orthW_RNA) {
    denom1 <- Pi_RNA * (W %*% t(W) %*% X_RNA %*% t(H_RNA) + L1 + L2*W)
} else {
    denom1 <- Pi_RNA * ((WH^(Beta-1) %*% t(H_RNA)) + L1 + L2*W)
}
```

Replace with weighted blend:

```r
denom1_std  <- WH^(Beta-1) %*% t(H_RNA)
denom1_orth <- W %*% t(W) %*% X_RNA %*% t(H_RNA)
denom1 <- Pi_RNA * ((1 - lambda_orthW) * denom1_std + lambda_orthW * denom1_orth + L1 + L2*W)
```

Apply identically in:
- `R/updateW_RNA2.R` Matrix mode
- `R/updateW_RNA2.R` List mode (per-chrom in the `Reduce` lambda)
- `R/updateW_RNA.R` (asymmetric Machima counterpart) — same change

## Tasks

1. **`R/Machima2.R` and `R/Machima.R`** — add `lambda_orthW = 0` to argument list (place near `orthW_RNA`). Pass through `.checkMachima2`/`.checkMachima` and `.initMachima2`/`.initMachima`. Pass through to `.updateW_RNA2`/`.updateW_RNA`.

2. **`R/checkMachima2.R` and `R/checkMachima.R`** — validate `is.numeric(lambda_orthW) && length == 1 && 0 <= lambda_orthW <= 1`. If `orthW_RNA = TRUE`, override `lambda_orthW = 1` and emit `warning("orthW_RNA=TRUE is deprecated; setting lambda_orthW=1. Use lambda_orthW directly.")`.

3. **`R/updateW_RNA2.R` and `R/updateW_RNA.R`** — change signatures to accept `lambda_orthW`, replace branch with blend (above). Drop the `orthW_RNA` parameter from these internal helpers (it's resolved at the public-API level).

4. **`man/Machima.Rd` and `man/Machima2.Rd`** — add `@param lambda_orthW Strength of W_RNA column orthogonality regularization. 0 = standard NMF, 1 = full orthogonal NMF, intermediate = weighted blend. (Default: 0)`.

5. **Tests** — add `tests/testthat/test-Machima2-lambda-orthW.R`:
   - `lambda_orthW = 0` reproduces `orthW_RNA = FALSE` exactly
   - `lambda_orthW = 1` reproduces `orthW_RNA = TRUE` exactly
   - `orthW_RNA = TRUE` triggers deprecation warning
   - intermediate values produce valid (non-NaN, non-negative) `W_RNA`
   - argument validation: `expect_error(Machima2(..., lambda_orthW = -0.1))`

6. **NEWS** — add `v1.4.0` entry to `inst/NEWS`:
   ```
   VERSION 1.4.0
   ------------------------
      o Added lambda_orthW in [0, 1] to Machima() and Machima2() for
        continuous control of W_RNA column orthogonality (replacing
        binary orthW_RNA toggle). Default 0 reproduces pre-1.4.0
        non-orthogonal behaviour.
      o orthW_RNA is now deprecated (still functional, emits warning).
   ```

7. **`DESCRIPTION`** — bump `Version: 1.4.0`.

## Out of scope

- Do **not** modify `orthH_RNA`, `orthT`, `orthH_Sym`. They retain their binary toggle. (Can be extended later if empirically motivated.)
- Do **not** add cross-modality orthogonality regularizers.

## Acceptance

- `R CMD check` clean.
- All existing tests pass with `lambda_orthW = 0` (default) — bit-exact identical to v1.3.0.
- New `test-Machima2-lambda-orthW.R` passes.
- Deprecation warning fires when `orthW_RNA = TRUE` is supplied.

Commit message: short subject describing parameter, body explaining motivation (binary toggle insufficient for tuning W column balance vs reconstruction quality, observed empirically in AIRE-experiments brain 100kb evaluation).
