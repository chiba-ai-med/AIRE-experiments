# Prompt: Implement `T_regularization = "low_rank"` for Machima2

## Audience

This prompt is for another Claude Code session that will modify the
**Machima R package** (https://github.com/kokitsuyuzaki/Machima),
v1.0.0 / SHA `c1d1804`. The user is the package maintainer.

## Context

`Machima2()` already supports two regularization strategies for the
free dense T matrix:

| `T_regularization` | implementation                                              |
|--------------------|-------------------------------------------------------------|
| `"none"`           | `T` is updated by `.updateT2` with no extra constraint      |
| `"frobenius_unit"` | per-iter rescale: `T[k] <- T[k]/c`, `H_Sym <- H_Sym * c²`   |
| `"l2"`             | `effective_L2_T = L2_T + lambda_T` is fed to `.updateT2`    |

Both attempt to prevent the `H_Sym` collapse symptom documented in
`notes/machima2_dense_t_regularization_prompt.md` — when `T` has
`l_k × n_k` (~ millions per chrom at 100kb) free parameters and `H_Sym`
has only `J²` (~ 49), the optimizer routes all magnitude into `T`.

The `frobenius_unit` path is a per-iter rescale heuristic; the `l2`
path requires `lambda_T` tuning. Option C from the original prompt —
**low-rank parametrization of T** — is a more principled fix: it
*structurally* bounds T's capacity, so no per-iter trick or hyperparameter
sweep is needed.

## Goal

Add `T_regularization = "low_rank"` as a third option, parametrizing

```
T[k] = U[k] · V[k]ᵀ                           # T[k]: l_k × n_k
       U[k]: l_k × r,  V[k]: n_k × r        # r = T_rank, user-supplied
```

Capacity per chrom is `(l_k + n_k) · r`, vs `l_k · n_k` in the dense
case. With `r << min(l_k, n_k)` (e.g. `r = 2J` or `r = J`), T no longer
has the slack to absorb all the magnitude away from H_Sym.

## API surface

Add one new argument to `Machima2()`:

```r
Machima2(
  ...,
  T_regularization = c("none", "frobenius_unit", "l2", "low_rank"),
  lambda_T = 0,
  T_rank   = NULL,         # NEW; required when T_regularization == "low_rank"
  ...
)
```

Defaults stay back-compat: `T_regularization = "none"`, `T_rank = NULL`.

Validation in `.checkMachima2`:

- When `T_regularization == "low_rank"`: `T_rank` must be a positive
  integer ≤ `min(l_k, n_k)` for every chrom (or for the matrix in
  matrix mode). Error message names the offending chrom and observed
  vs allowed bounds.
- When `T_regularization != "low_rank"`: `T_rank` must be `NULL` or
  ignored with a `warning()` ("`T_rank` ignored when
  `T_regularization` is not 'low_rank'").
- `T_regularization == "low_rank"` is incompatible with `fixT = TRUE`
  (because `fixT` skips the T update entirely, including the low-rank
  reparam). Error if both set.
- `T_regularization == "low_rank"` is also incompatible with the
  user passing a non-NULL `T` argument that is full-rank dense — the
  caller would need to factor it themselves or accept that it's
  re-initialized. Either error out, or accept the user's `T` and
  initialize `U`, `V` from a rank-`r` truncated SVD of it.

## Implementation sketch

### Initialization (`R/initMachima2.R`)

Add a branch at the end of `.initMachima2_List` and `.initMachima2_Matrix`:

```r
if(T_regularization == "low_rank"){
    if(is.null(int$T)){
        # Random-init U, V (non-negative)
        U <- lapply(seq_along(X_RNA), function(k){
            matrix(runif(nrow(X_Epi[[k]]) * T_rank), ncol = T_rank)
        })
        V <- lapply(seq_along(X_RNA), function(k){
            matrix(runif(ncol(X_RNA[[k]]) * T_rank), ncol = T_rank)
            # actually V[k] should be n_k × r, where n_k = nrow(X_RNA[[k]])
        })
    }else{
        # User supplied a dense T -- factor it
        U <- lapply(int$T, function(t){
            sv <- svd(t, nu = T_rank, nv = T_rank)
            sv$u * sqrt(sv$d[seq_len(T_rank)])[col(sv$u)]
        })
        V <- ...   # symmetric
    }
    # T[k] is the materialized product, kept consistent with U,V at all times
    int$T <- lapply(seq_along(U), function(k) U[[k]] %*% t(V[[k]]))
    int$U <- U
    int$V <- V
}
```

Carry `U` and `V` through the iteration loop in addition to (or instead
of) the materialized `T`. Easiest path: keep `T = U Vᵀ` materialized
each iter for the rest of the code (which references `T` directly), but
update `U` and `V` instead of `T`.

### Update step (`R/updateT2.R` or new `R/updateT2_lowrank.R`)

Replace the dense `.updateT2` call with two multiplicative updates:

```r
# Multiplicative NMF-style update for U[k], V[k] given T[k] = U[k] V[k]ᵀ
# Loss: ‖X_Epi[k] - (U[k] V[k]ᵀ W H_Sym Wᵀ V[k] U[k]ᵀ)‖_F²
# Derive ∂loss/∂U[k] and ∂loss/∂V[k] and apply standard
# multiplicative-update form: U <- U * (numerator / denominator).
```

The derivation is mechanical from the dense case — `U[k]` plays the role
of the left factor, `V[k] W H_Sym Wᵀ V[k]ᵀ` plays the role of the
"effective core". Implement with care for numerical stability
(`pmax(denominator, eps)`).

After updating `U` and `V`, materialize `T[k] <- U[k] %*% t(V[k])` so
the rest of the loop (which uses `T` as the model variable) works
unchanged.

### Storage on the result list

`Machima2()` returns `T` as part of the result; for backward compat,
keep returning the materialized `T = U Vᵀ`. Optionally also return
`U` and `V` as `res$T_factors$U`, `res$T_factors$V` so users can
inspect or warm-start.

## Acceptance test

1. **API**: `Machima2(..., T_regularization = "low_rank", T_rank = 14)`
   on a list-mode input (J=7, K=22 chroms) runs to completion. With
   `T_rank = 7`, ditto.
2. **Rank**: After fit, `qr(res$T[[k]])$rank` ≤ `T_rank` for every k.
3. **H_Sym informative**: On AIRE-experiments / brain 100kb joint+Tdense
   with `T_regularization = "low_rank"`, `T_rank = 14`,
   `diag(H_Sym)` values are O(1) (i.e. the H_Sym collapse symptom is
   gone, matching the `frobenius_unit` acceptance criterion).
4. **Hi-C reconstruction not regressed**: `rel_frob` stays in 0.12–0.20
   (matches the dense-T baseline; should not be much worse since
   `r = 14` is still well above the J=7 informative rank).
5. **fixT compat**: `Machima2(..., fixT = TRUE, T_regularization = "low_rank")`
   throws a clear error.
6. **`T_rank = NULL` validation**: `Machima2(..., T_regularization = "low_rank")`
   without `T_rank` throws a clear error naming the missing arg.
7. **Existing tests unchanged**: All current tests
   (`test-Machima2-init.R`, `test-Machima2-Treg.R`, `test-dimnames-na.R`)
   pass unchanged.

## Out of scope

- Optimizing the multiplicative update (e.g., closed-form alternating
  least squares for U/V) — a multiplicative-update form is fine for
  v1.0.x.
- Auto-selecting `T_rank` (e.g. via reconstruction-error elbow). User
  supplies it explicitly.
- Adding low-rank to identity T — meaningless (identity T is already
  rank `n_k`).
- Removing the existing `"frobenius_unit"` and `"l2"` modes.

## Useful pointers

- Original Option C spec: `notes/machima2_dense_t_regularization_prompt.md`
  (the section beginning `### Option C (most invasive): low-rank T`)
- Matching `frobenius_unit` list-mode bug being addressed in parallel:
  `notes/machima2_frobenius_unit_listmode_fix_prompt.md`
- Diagnostic that motivates the regularization need:
  `output/diagnostic_brain_100000_hsym_summary.csv`
