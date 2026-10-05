# Machima v1.3.0 — Expose `lambda_balance` for explicit X_RNA / X_Epi loss balancing

You are working on the upstream R package **`kokitsuyuzaki/Machima`**. The package provides `Machima2()` (and `Machima()`), joint NMTF for cell-type deconvolution. Two losses are summed:

```
L = L_RNA(W, H_RNA, X_RNA)  +  L_Epi(W, H_Sym, T, X_Epi)
```

## What already exists

`R/Machima-internal.R` defines

```r
.weight <- function(X){ 1 / sum(X^2) }
```

and `.initMachima` / `.initMachima2` automatically compute

```r
Pi_RNA <- .weight(X_RNA)        # = 1 / ||X_RNA||²_F
Pi_Epi <- .weight(X_Epi)        # = 1 / ||X_Epi||²_F
```

These scalars are passed through every `.update*` helper and multiply both numerator and denominator of the MU update. So **Frobenius normalisation is already automatic and the two terms are weighted equally** by default. Good — but there is no user-facing knob to deviate from "equal".

## The motivation

In paired scATAC + Hi-C deconvolution, even with Frobenius-normalised loss, the shared `W_RNA` is in practice optimised more strongly by the X_RNA term because:

- Σ X_RNA gradient mass: O(n × m) — many cells, dense
- Σ X_Epi gradient mass: O(l²) — single bulk Hi-C, rich row-column structure but only one matrix per chromosome

Empirically we observed that joint factorisation collapses `H_Sym` toward a rank-1 outer product `d · dᵀ` (one shared contact pattern weighted differently per component) — the model finds it cheap to satisfy X_Epi by absorbing magnitudes into a few entries of `H_Sym` while leaving `W_RNA` slaved to RNA. Letting users **tilt the balance toward Epi** (or RNA) on demand is a clean lever to investigate this.

## The proposed parameter

Add **`lambda_balance ∈ [0, 1]`** to both `Machima()` and `Machima2()`:

- `lambda_balance = 0`   → 100% RNA loss, X_Epi ignored (Pi_Epi = 0 effectively)
- `lambda_balance = 1`   → 100% Epi loss, X_RNA ignored
- `lambda_balance = 0.5` → equal, **identical to current default**

**Default value: `0.5`** — preserves backward compatibility exactly.

## Implementation plan

### 1. Function signatures

In `R/Machima.R` and `R/Machima2.R`, add `lambda_balance = 0.5` to the argument list (place near `Pi_RNA` / similar weighting concepts; suggested position: just before `T_regularization`).

### 2. Argument check (`R/checkMachima.R` and `R/checkMachima2.R`)

Add:

```r
stopifnot(is.numeric(lambda_balance), length(lambda_balance) == 1,
          lambda_balance >= 0, lambda_balance <= 1)
```

Pass `lambda_balance` through from the public function to `.initMachima*`.

### 3. `.initMachima*` weight construction

Replace the current

```r
Pi_RNA <- .weight(X_RNA)
Pi_Epi <- .weight(X_Epi)
```

with

```r
Pi_RNA <- 2 * (1 - lambda_balance) * .weight(X_RNA)
Pi_Epi <- 2 *      lambda_balance  * .weight(X_Epi)
```

(The factor of 2 makes `lambda_balance = 0.5` yield exactly the v1.2.0 weights — `2 * 0.5 = 1`, so `Pi_RNA = .weight(X_RNA)`, identical to before.)

For list mode (matrix list), apply the same scalar multiplier to the per-element list:

```r
Pi_RNA <- lapply(X_RNA, function(x) 2 * (1 - lambda_balance) * .weight(x))
Pi_Epi <- lapply(X_Epi, function(x) 2 *      lambda_balance  * .weight(x))
```

### 4. Edge cases

- `lambda_balance = 0`: `Pi_Epi = 0`. Verify update rules don't divide by zero. Inspect `.updateW_RNA2` denom2 (which involves `Pi_Epi * (...)`) — when `Pi_Epi = 0`, denom2 = 0, and num2 = 0; their contributions to `(num1 + num2) / (denom1 + denom2)` become 0/0 only if num1 + num2 = 0 too. Add a `pseudocount`-equivalent guard: e.g., clamp `Pi_RNA` and `Pi_Epi` to `>= .Machine$double.eps` if either is exactly 0. Or document that exact endpoints `lambda_balance ∈ {0, 1}` are unsupported and require `1e-6` margin. Document whichever you choose in `@param`.
- Default-only callers (no `lambda_balance` argument supplied) get exactly v1.2.0 behaviour by construction. No regression.

### 5. roxygen documentation

In `R/Machima2.R`:

```r
#' @param lambda_balance Balance between X_RNA and X_Epi loss terms in the
#'   joint objective. Both terms are first Frobenius-normalised
#'   (`Pi_RNA = 1/||X_RNA||²`, `Pi_Epi = 1/||X_Epi||²`); `lambda_balance`
#'   then tilts between them: 0 = X_RNA only, 1 = X_Epi only,
#'   0.5 = equal (default; matches pre-1.3.0 behaviour). Useful when one
#'   modality dominates the joint optimisation. (Default: 0.5)
```

Identical wording in `R/Machima.R` (referring to `H_Epi` instead of `H_Sym` if relevant in the surrounding docs).

Add a sentence under `@details`:

> "By default both modality losses are Frobenius-normalised and contribute equally. Increasing `lambda_balance` toward 1 reweights the joint optimisation toward fitting `X_Epi` (epigenome) more strictly — useful when investigating cell-type-specific contact structure that can be drowned out by the larger `X_RNA` gradient mass."

### 6. Tests

Add `tests/testthat/test-Machima2-loss-balance.R`:

- **Default equivalence**: `Machima2(X_RNA, X_Epi, ..., lambda_balance = 0.5)` produces output **identical** to `Machima2(X_RNA, X_Epi, ...)` (no argument). Use `expect_equal(res_a$W_RNA, res_b$W_RNA)` etc. with same RNG seed.
- **Endpoint behaviour at `lambda_balance ≈ 1`**: with `lambda_balance = 0.99`, the X_Epi reconstruction error (computed per `.recErrors2`) should be lower than with `lambda_balance = 0.01`, on a fixed synthetic problem where ground-truth signal exists in X_Epi.
- **Endpoint behaviour at `lambda_balance ≈ 0`**: symmetric — X_RNA reconstruction better.
- **Argument validation**: `expect_error(Machima2(..., lambda_balance = -0.1))`, `expect_error(Machima2(..., lambda_balance = 1.5))`.

### 7. NEWS

Add a v1.3.0 entry to `inst/NEWS`:

```
VERSION 1.3.0
------------------------
   o Added lambda_balance ∈ [0, 1] to Machima() and Machima2() to control
     the relative weight of X_RNA vs X_Epi loss terms in the joint
     objective. Default 0.5 is identical to pre-1.3.0 behaviour
     (Frobenius-normalised equal weighting via existing Pi_RNA/Pi_Epi).
```

(If the diagonal-MU change ships in the same v1.3.0, list both bullets together.)

### 8. DESCRIPTION

Bump `Version: 1.3.0` if not already done.

## Out of scope

- Do **not** remove `.weight()` or change its semantics — it is the per-modality Frobenius normaliser and stays intact.
- Do **not** add other loss-balance schemes (e.g. KL-divergence-based normalisation, learnable weights). Just the single scalar `lambda_balance`.
- Do **not** add per-chromosome balancing parameters in list mode. The lambda is global.

## Acceptance criteria

- `R CMD check` clean.
- All existing tests pass unchanged (default `lambda_balance = 0.5` reproduces v1.2.0 behaviour).
- New `test-Machima2-loss-balance.R` passes.
- `inst/NEWS` and `DESCRIPTION` reflect v1.3.0.
- The change is exposed identically in `Machima()` (the asymmetric variant) and `Machima2()` (the symmetric variant).

Commit message convention: short subject (≤72 chars) describing the parameter, body explaining the existing `.weight()` mechanism and why a user-facing knob is useful for diagnosing optimisation balance issues.
