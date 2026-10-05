# Machima v1.3.0 — Replace ad hoc projection with proper constrained MU for `H_Sym_structure="diagonal"`

You are working on the upstream R package **`kokitsuyuzaki/Machima`** (GitHub-only). The package provides `Machima2()`, a symmetric NMTF for cell-type deconvolution of bulk Hi-C using paired scATAC + scRNA multiome as the cell-type basis. Model:

- `X_RNA  ≈ W_RNA · H_RNA`  (per cell)
- `X_Epi  ≈ G · H_Sym · Gᵀ`,  where `G = T · W_RNA` (per chromosome; `H_Sym` shared)

## Background

In **v1.2.0** we added the option `H_Sym_structure = c("symmetric", "diagonal")`. Diagonal mode (sym-CP / sym-PARAFAC) was supposed to make the deconvolution principled: bulk Hi-C contacts only occur within single cells of one type, so cross-component (off-diagonal) `H_Sym[i,j]` entries have no biological meaning. The diagonal model is

```
X_Epi[k]  ≈  Σ_i  h_i · g_i^{(k)} · g_i^{(k)ᵀ},   h ∈ ℝ^J_+
```

where `g_i^{(k)} = (T_k · W_RNA)[:,i]` is the i-th component vector on chromosome k, and `H_Sym = diag(h)`.

## The problem

Empirical evaluation on brain 100 kb (joint × Tidentity) showed diagonal mode collapses fitting catastrophically:

| metric | symmetric (sym-Tucker) | diagonal (sym-CP) |
|---|---|---|
| RecError | 3.64 | 13.42 |
| mean_rel_frob_per_chr | 0.123 | **0.995** (≈ random) |
| OVERALL pearson (matched) | 0.508 | 0.056 |

Some loss is expected (Tucker ≥ CP in fit), but **8× degradation is too large**. Investigation of the v1.2.0 implementation reveals the cause: diagonal mode is implemented as **ad hoc post-MU projection**:

```r
# R/Machima2.R l.166, l.184  (both horizontal and tri-factorisation modes)
H_Sym <- .updateH_Sym(X_Epi, W_RNA, H_Sym, T, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
if (H_Sym_structure == "diagonal") H_Sym <- diag(diag(H_Sym))
```

This computes the **full J×J MU update** (`.updateH_Sym` operates on the unconstrained symmetric H_Sym) and then zeros off-diagonals. The descent direction is **for the unconstrained problem**, not the diagonal-constrained problem. After projection there is no convergence guarantee, and per-iteration progress is far below what proper constrained MU achieves.

**Fix**: parameterize `h ∈ ℝ^J_+` directly and derive a proper MU rule on `h`.

## Math (constrained MU on the diagonal manifold)

For β-divergence loss

```
L(h) = Σ_k  d_β( X_Epi^{(k)}, R^{(k)}(h) ),   R^{(k)}(h) = Σ_i h_i g_i^{(k)} g_i^{(k)ᵀ}
```

the gradient is

```
∂L/∂h_i = Σ_k g_i^{(k)ᵀ} ( R^{(k)}^{β-1}  −  X_Epi^{(k)} ⊙ R^{(k)}^{β-2} ) g_i^{(k)}
```

(elementwise power; `⊙` is Hadamard.) Multiplicative split into positive numerator and denominator:

```
num_i   = Σ_k  g_i^{(k)ᵀ} ( X_Epi^{(k)} ⊙ R^{(k)}^{β-2} ) g_i^{(k)}     # from -∂L/∂h_i
denom_i = Σ_k  g_i^{(k)ᵀ}              R^{(k)}^{β-1}     g_i^{(k)}     # from +∂L/∂h_i
```

MU rule with relaxation exponent `ρ(β)` (matching existing convention in `R/updateH_Sym.R`):

```
h_i  ←  h_i · (num_i / denom_i)^{ρ(β)}
```

L1/L2 regularisation translates as in standard NMF:

- L1: `denom_i ← denom_i + L1_H_Sym`
- L2: `denom_i ← denom_i + L2_H_Sym · h_i`

`Pi_Epi` (per-modality weight) multiplies both `num` and `denom` element-wise (or in list mode, per-k). `orderReg` is meaningful (orders h descending across i) and should be honoured if currently honoured for the unconstrained case. `orthH_Sym` is **trivially satisfied** for diagonal H_Sym and should be silently ignored (or warned at most once at init).

For β = 2 (Frobenius) the formula simplifies to:

```
num_i   = Σ_k  g_i^{(k)ᵀ} X_Epi^{(k)} g_i^{(k)}
denom_i = Σ_k  Σ_j  h_j · ⟨g_i^{(k)}, g_j^{(k)}⟩²
```

The second form is O(J²) per chromosome and avoids forming `R^{(k)}` explicitly — exploit this for speed.

## Tasks

1. **Add helper `R/updateH_Sym_diag.R`** with two functions (mirroring the existing `.updateH_Sym` / `.updateH_Sym_HZL` pair in `R/updateH_Sym.R`):
   - `.updateH_Sym_diag(X_Epi, W_RNA, h, T, J, Beta, L1, L2, orderReg, root, Pi_Epi)` — tri-factorisation mode
   - `.updateH_Sym_diag_HZL(X_GAM, W_RNA, h, J, Beta, L1, L2, orderReg, root, Pi_Epi)` — horizontal mode (use the same `X_GAM` / `H_Sym ≈ aux` shortcut already used by `.updateH_Sym_HZL`)

   Both accept a length-J numeric vector `h` and return an updated length-J numeric vector. Internally they may use the J×J matrix form `diag(h)` for matrix algebra but must NOT update off-diagonals.

   For both Matrix and List modes: accumulate `num` and `denom` over `k` (for List input) or compute once (for Matrix input). Keep code paths consistent with the existing `.updateH_Sym` style in `R/updateH_Sym.R` (look at how that function dispatches Matrix vs List).

2. **Modify `R/Machima2.R`** to dispatch on `H_Sym_structure`:

   Replace
   ```r
   if (!fixH_Sym) {
       H_Sym <- .updateH_Sym(...)
       if (H_Sym_structure == "diagonal") H_Sym <- diag(diag(H_Sym))
   }
   ```
   with
   ```r
   if (!fixH_Sym) {
       if (H_Sym_structure == "diagonal") {
           h <- .updateH_Sym_diag(X_Epi, W_RNA, diag(H_Sym), T, J, Beta,
               L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi)
           H_Sym <- diag(h)
           # preserve dimnames if previously assigned
       } else {
           H_Sym <- .updateH_Sym(X_Epi, W_RNA, H_Sym, T, J, Beta,
               L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
       }
   }
   ```
   Apply the same change in both `horizontal` and tri-factorisation branches of the iteration loop.

3. **Verify `.initMachima2` still does its diagonal projection at init time** (lines ~93 and ~213 in `R/initMachima2.R`). Init-time projection is correct (it just sets the starting point on the diagonal manifold). Leave that alone.

4. **Update `man/Machima2.Rd`** `@param H_Sym_structure` description to note that "diagonal" now uses a constrained MU (rather than ad hoc post-projection). Add a sentence under `@details`:

   > "When `H_Sym_structure = "diagonal"`, the diagonal entries `h ∈ ℝ^J_+` are updated by a proper constrained MU rule derived from the J-parameter loss, not by post-MU projection."

5. **Tests** — add `tests/testthat/test-Machima2-Hstruct-mu.R`:

   a. **Diagonal preservation across iterations**: run with `H_Sym_structure="diagonal"`, verify `all(abs(res$H_Sym[upper.tri(res$H_Sym)]) < 1e-12)` after every iteration (not only at the end). One way: set `num.iter=10`, then run 10 separate 1-iter calls with `init_H_Sym = previous H_Sym`, asserting at each step.

   b. **Monotonic RecError decrease (within Frobenius)**: with `Beta=2`, `RecError` must be non-increasing across iterations (modulo floating-point tolerance ~1e-10). The ad hoc projection does NOT satisfy this guarantee; the proper MU does.

   c. **Numerical equivalence with toy rank-1 case**: when J=1, diagonal and symmetric updates are identical. Test that both produce the same output.

   d. **Fitting on a synthetic sym-CP signal**: generate `G ∈ ℝ^{15×3}_+` and `h ∈ ℝ^3_+`, set `X_Epi = G diag(h) Gᵀ + small noise`, run with `H_Sym_structure="diagonal"`. Assert recovered `H_Sym` is close to ground-truth `diag(h)` (matched up to permutation of components).

6. **NEWS** — add a v1.3.0 entry to `inst/NEWS`:

   ```
   VERSION 1.3.0
   ------------------------
      o H_Sym_structure="diagonal" now uses a proper constrained
        multiplicative update (acting directly on the J diagonal entries
        h ∈ R^J), rather than ad hoc post-MU projection of a full J×J
        update. Yields large fit improvements on Hi-C data; previous
        ad hoc behaviour is no longer reachable.
      o Added .updateH_Sym_diag / .updateH_Sym_diag_HZL helpers.
   ```

7. **DESCRIPTION** — bump `Version: 1.3.0`.

## Out of scope

- Do **not** change the symmetric (sym-Tucker) MU rule. It is correct as-is.
- Do **not** alter `H_Sym_structure` defaults. Default remains `"symmetric"` (it fits better on real Hi-C; users opt into diagonal).
- The `init_H_Sym` warning logic in `.initMachima2` (warns if user-supplied init has non-zero off-diagonals when `H_Sym_structure="diagonal"`) stays.

## Acceptance criteria

- All existing `tests/testthat/*.R` still pass (run `R CMD check` clean — no new ERROR/WARNING).
- New `test-Machima2-Hstruct-mu.R` passes.
- `inst/NEWS` and `DESCRIPTION` updated to v1.3.0.
- No new external dependencies.

Commit message convention: short subject (≤72 chars) describing the substantive change, followed by a body explaining what was wrong with the v1.2.0 ad hoc projection and how the proper constrained MU resolves it.
