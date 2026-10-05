# Task: Add `init_*` parameters to `Machima2()` for two-stage workflows

## Repository
https://github.com/kokitsuyuzaki/Machima  (R package, GitHub-only — not on CRAN/Bioconductor)

## Background

`Machima2()` is a multi-modal NMF that jointly factorizes paired
single-cell counts (`X_RNA`, n × m) and bulk Hi-C contact maps
(`X_Epi`, l × l symmetric) under the model

```
X_RNA ≈ W_RNA · H_RNA
X_Epi ≈ (T · W_RNA) · H_Sym · (T · W_RNA)ᵀ
```

with `T` either fixed (identity) or learned (dense). The shared factor
`W_RNA` couples the two modalities.

In `R/Machima2.R` the function already exposes `fixW_RNA`, `fixH_RNA`,
`fixT`, `fixH_Sym` flags but **does NOT accept user-supplied initial
values for `W_RNA`, `H_RNA`, `H_Sym`** (only `T` can be passed in).
Initialization is controlled by an enum
`init = c("Random", "RandomEpi", "RandomRNA")` and the resulting init
matrices are constructed inside `.initMachima2_List` / `.initMachima2_Matrix`
via `.returnBestNMF()` which calls `nnTensor::NMF()` with hard-coded
`num.iter = 30, algorithm = "Frobenius"`.

## Goal

Add user-supplied initialization for **all** factors so callers can do
two-stage / staged-fitting workflows like:

1. Run `nnTensor::NMF` themselves (with their own `algorithm`,
   `num.iter`, restarts) on `X_RNA` to obtain `W_RNA_0`, `H_RNA_0`.
2. Pass `W_RNA_0` and `H_RNA_0` into `Machima2()` as init, set
   `fixW_RNA = TRUE, fixH_RNA = TRUE`, and let Machima2 only fit
   `H_Sym` (and `T` if dense) against Hi-C.

This eliminates the need for monkey-patches like
`assignInNamespace(".returnBestNMF", ...)` that downstream users have
been forced to use.

## API change

Add four new named arguments to `Machima2()` (defaulting to `NULL` for
backward compatibility — current callers should be unaffected):

```r
Machima2(
  ...,
  init_W_RNA = NULL,  # list of (n_k × J) matrices, OR an (n × J) matrix
                      # — must match the shape contract of X_RNA
  init_H_RNA = NULL,  # (J × m) matrix
  init_H_Sym = NULL,  # (J × J) symmetric matrix
  # init_T already exists as the `T` argument — keep that name unchanged
  ...
)
```

Behaviour:

- If `init_W_RNA` is supplied, skip the random/RandomEpi/RandomRNA init
  for `W_RNA` and use the user-supplied value as the starting point.
- Same for `init_H_RNA` and `init_H_Sym`.
- The existing `init` enum still controls how the **un-supplied** factors
  are initialized (so `init_W_RNA = my_W, init = "RandomRNA"` would use
  `my_W` for W and random for the rest).
- Combines naturally with `fixW_RNA = TRUE` etc: if both are set, the
  factor is initialized from the user's matrix and never updated.
- When the input is a list (per-chromosome mode), `init_W_RNA` must be
  a list of the same length and per-element row counts must match
  `nrow(X_RNA[[k]])`.

## Validation

Extend `.checkMachima2()` to validate the new arguments:

- Type check: NULL, matrix, or list-of-matrices as appropriate.
- Dimension match against `X_RNA` / `X_Epi` / `J`.
- For `init_H_Sym`: must be J×J and symmetric (`isSymmetric` after
  `as.matrix()`).
- Clear error messages naming the offending argument and observed vs
  expected dimensions.

## Wiring inside `.initMachima2_List` / `.initMachima2_Matrix`

After the existing `init`-based construction of `W_RNA`, `H_RNA`,
`H_Sym`, overwrite each one with the user-supplied value if non-NULL.
Do this *after* the random-init branch so the random matrices act as
fallback for any factor the user did not specify. Do NOT recompute
`Pi_RNA` / `Pi_Epi` from the supplied factors — they're derived from
the data only.

## Secondary fix: `.returnBestNMF` restart bug

The current implementation:

```r
.returnBestNMF <- function(X, J) {
    outs <- lapply(seq(1), function(x) {        # ← seq(1) = c(1), 1 restart
        NMF(X, J = J, num.iter = 30, algorithm = "Frobenius")
    })
    bestfit <- unlist(lapply(outs, function(x) rev(x$RecError)[1]))
    bestfit <- which(bestfit == min(bestfit))[1]
    outs[[bestfit]]
}
```

does only 1 restart despite the `lapply` + `bestfit` selection
suggesting multiple. Fix by accepting a `n_restart` (and `num_iter`,
`algorithm`) argument, defaulting to current behaviour for back-compat:

```r
.returnBestNMF <- function(X, J, n_restart = 1L,
                           num_iter = 30L,
                           algorithm = "Frobenius") {
  outs <- lapply(seq_len(n_restart), function(x) {
    nnTensor::NMF(X, J = J, num.iter = num_iter, algorithm = algorithm)
  })
  ...
}
```

Plumb these through to `Machima2()` as new args
`nmf_init_n_restart`, `nmf_init_num_iter`, `nmf_init_algorithm` so
users can control the internal `RandomEpi`/`RandomRNA` NMF runs without
having to pre-compute and pass `init_W_RNA`/`init_H_RNA` themselves.

## Documentation

Update the Roxygen `@param` block on `Machima2()` for all new
arguments. In `@examples`, add one example showing the two-stage
pattern:

```r
\dontrun{
# Stage A: pure RNA NMF with user-controlled algorithm + restarts
nmf_res <- nnTensor::NMF(do.call(rbind, X_RNA),
                         J = 5, num.iter = 200,
                         algorithm = "KL")
W0 <- # split nmf_res$U back into per-chrom list (helper .Mat2List)
H0 <- nmf_res$V

# Stage B: transfer into Machima2, freeze the RNA factors
res <- Machima2(
  X_RNA = X_RNA, X_Epi = X_Epi,
  init_W_RNA = W0, init_H_RNA = H0,
  fixW_RNA = TRUE, fixH_RNA = TRUE,
  J = 5, num.iter = 100
)
}
```

If `.Mat2List` is internal-only (`:::`), expose it as a user-facing
helper (`Mat2List`) since callers will need it for the per-chrom split.

## Acceptance criteria

1. Passing `init_W_RNA = res_old$W_RNA` from a previous fit reproduces
   the same starting iteration (verify by setting `num.iter = 1` and
   inspecting `RecError[1]`).
2. With all four init_* supplied AND all four fix* set, the algorithm
   does no updates and `res$W_RNA` etc. equal the supplied matrices
   exactly.
3. Existing test suite passes unchanged (back-compat for current
   call signature).
4. New test: two-stage workflow (Stage A `nnTensor::NMF` then Stage B
   `Machima2(..., init_W_RNA=..., fixW_RNA=TRUE, fixH_RNA=TRUE)`)
   produces a `res$W_RNA` element-wise identical to `init_W_RNA`.
5. `R CMD check` clean, no new NOTEs.

## Out of scope

- Changing the update rules themselves.
- Adding per-modality `Beta` (currently `Beta` is global to both
  X_RNA and X_Epi reconstruction — that's a separate, larger refactor).
- Performance optimization.

## Files likely touched

- `R/Machima2.R` — function signature + init wiring
- `R/initMachima2.R` (or wherever `.initMachima2_List` /
  `.initMachima2_Matrix` live) — accept and propagate init_*
- `R/checkMachima2.R` — validate init_*
- `R/returnBestNMF.R` — restart/algorithm/iter args
- `R/Mat2List.R` (new export) — if not already user-facing
- `man/Machima2.Rd` — regenerated by Roxygen
- `tests/testthat/test-Machima2-init.R` — new tests
- `NEWS.md` — note the new args
