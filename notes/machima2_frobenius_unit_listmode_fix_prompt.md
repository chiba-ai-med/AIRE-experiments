# Prompt: Fix `frobenius_unit` rescale to be math-correct in list mode (per-chrom T)

## Audience

This prompt is for another Claude Code session that will modify the
**Machima R package** (https://github.com/kokitsuyuzaki/Machima),
v1.0.0 / SHA `c1d1804`. The user is the package maintainer.

## Context

Machima2's symmetric model is

```
X_RNA[k] ≈  W_RNA[k] · H_RNA
X_Epi[k] ≈ (T[k] · W_RNA[k]) · H_Sym · (T[k] · W_RNA[k])ᵀ        # H_Sym shared across all k
```

In commit `c1d1804`, `T_regularization = "frobenius_unit"` was added to
prevent H_Sym from collapsing to ~1e-12 when T is a free dense matrix.
The current implementation, in `R/Machima2-internal.R`:

```r
.frobNormT <- function(T){
    if(is.matrix(T)){
        nrm <- norm(T, "F")
        list(norms = nrm, mean_norm_sq = nrm^2)
    }else{
        norms <- sapply(T, function(t) norm(t, "F"))
        list(norms = norms, mean_norm_sq = mean(norms)^2)
    }
}

.rescaleT <- function(T, frob){
    if(is.matrix(T)){
        T / frob$norms
    }else{
        lapply(seq_along(T), function(x) T[[x]] / frob$norms[x])   # ← per-chrom scalar
    }
}
```

and in `R/Machima2.R` main loop:

```r
if(T_regularization == "frobenius_unit"){
    frob  <- .frobNormT(T)
    T     <- .rescaleT(T, frob)
    H_Sym <- H_Sym * frob$mean_norm_sq                              # ← mean(c_k)^2
}
```

## Symptom (problem to solve)

When AIRE-experiments / brain 100 kb is rerun with
`T_regularization = "frobenius_unit"` for the joint+Tdense stage, the
fit crashes at iteration 13:

```
Error in while ((RelChange[iter] > thr) && (iter <= num.iter)) :
  missing value where TRUE/FALSE needed
```

(Source: `logs/machima/run_machima2_brain_joint_Tdense_100000.log` in
the AIRE-experiments repo.) `RelChange[iter]` is `NA`, which traces back
to `RecError[iter]` becoming `NaN`/`Inf`.

## Why this happens — math is broken in list mode

Each `T[k]` is divided by its **own** Frobenius norm `c_k = ‖T[k]‖_F`,
so the chrom-k reconstruction shrinks by `1/c_k²`:

```
(T_k / c_k)·W·H_Sym·Wᵀ·(T_k / c_k)ᵀ  =  (1/c_k²)·(T_k·W·H_Sym·Wᵀ·T_kᵀ)
```

To preserve the reconstruction *exactly*, H_Sym for chrom k would need
to be multiplied by `c_k²`. But `H_Sym` is **shared across all k**, so
any single global compensation can be exact for at most one chrom. The
current implementation uses `mean(c_k)²`, which is exact only when all
`c_k` are identical. For unequal `c_k` (the realistic case — chr1 vs
chr19 differ in size by an order of magnitude), the per-chrom
reconstruction error spikes immediately after each rescale. The next
multiplicative T update then propagates that spike, and within ~10–15
iterations the per-element products go to 0 / Inf and produce NaN.

Additional failure mode: if any `c_k` becomes 0 (e.g. a chrom whose `T[k]`
collapses to zero during a multiplicative update), `T[k] / 0 = NaN/Inf`
on the very next rescale step, regardless of compensation strategy.

## Required change

Two independent fixes — both should land. Either alone closes the NaN,
but only both together yield mathematically clean behavior.

### Fix A (mandatory): use a single global scalar in list mode

Change the list branch of `.frobNormT` and `.rescaleT` to compute **one
scalar `c`** and apply it to every `T[k]` (and to `H_Sym` as `c²`).
Recommended choice: RMS of per-chrom norms,

```
c = sqrt( mean_k c_k² ) = sqrt( sum_k ‖T[k]‖_F² / K )
```

Geometric mean `c = (∏_k c_k)^(1/K)` is also acceptable; arithmetic mean
is not (RMS keeps the *aggregate* Frobenius norm `sqrt(∑ ‖T[k]‖_F²)` at
the constant `sqrt(K)` after rescale, which is the natural per-chrom-1
target).

After rescale,

```
T[k]   <- T[k] / c     for all k
H_Sym  <- H_Sym * c²
```

This preserves `(T[k]·W) · H_Sym · (T[k]·W)ᵀ` exactly for every chrom.

### Fix B (mandatory): eps guard

Add a small floor to `c` to prevent division by zero:

```r
eps <- sqrt(.Machine$double.eps)   # ~ 1.49e-8
c   <- max(c, eps)
```

This must be applied to the global scalar from Fix A. (If you keep a
per-chrom branch as a fallback for any reason, the same guard applies
per element.)

## Suggested implementation

```r
.frobNormT <- function(T){
    eps <- sqrt(.Machine$double.eps)
    if(is.matrix(T)){
        nrm <- max(norm(T, "F"), eps)
        list(scalar = nrm, scalar_sq = nrm^2)
    }else{
        norms <- sapply(T, function(t) norm(t, "F"))
        # RMS over chroms; eps-guarded
        c <- sqrt(mean(norms^2))
        c <- max(c, eps)
        list(scalar = c, scalar_sq = c^2)
    }
}

.rescaleT <- function(T, frob){
    if(is.matrix(T)){
        T / frob$scalar
    }else{
        lapply(T, function(t) t / frob$scalar)
    }
}
```

Main loop becomes:

```r
if(T_regularization == "frobenius_unit"){
    frob  <- .frobNormT(T)
    T     <- .rescaleT(T, frob)
    H_Sym <- H_Sym * frob$scalar_sq
}
```

The list-element field name `norms` (plural) goes away — there's only
one scalar now in either branch. Rename to `scalar` / `scalar_sq` to
make the intent clear, or keep the `norms` field name for backward
compatibility if anything internal touches it (it shouldn't, since it's
package-internal `.`-prefixed).

## Acceptance test

The user will rerun AIRE-experiments / brain 100 kb with
`T_regularization = "frobenius_unit"` and check:

1. **No NaN crash.** All 6 Dense T stages (joint, transferFlog,
   transferKraw, supervisedHfix, supervisedGNMF, supervisedWinit)
   complete the configured 30 iterations without hitting the
   `missing value where TRUE/FALSE needed` error.
2. **Reconstruction equality.** Add a unit test that fits Machima2 with
   `T_regularization = "none"` and `T_regularization = "frobenius_unit"`
   on the same random seed, and confirms that the reconstruction
   `(T·W)·H_Sym·(T·W)ᵀ` is **identical to within numerical tolerance**
   immediately after the rescale step. (The optimization trajectory
   will differ, but the rescale itself must be model-preserving.)
3. **Aggregate norm stable.** After the rescale,
   `sqrt(sum_k ‖T[k]‖_F²) ≈ sqrt(K)` (within 1e-6) for list mode, and
   `‖T‖_F ≈ 1` for matrix mode.
4. **`H_Sym` magnitude.** For brain 100kb joint+Tdense,
   `diag(H_Sym)` should be O(1) — not O(1e-12) (the bug we're trying
   to fix in the first place) and not exploding (the new bug we're
   trying to avoid).
5. **Identity T unaffected.** `T_regularization` is silently ignored
   when `fixT = TRUE`. Existing identity-T tests pass unchanged.

## Out of scope

- Switching to per-chrom `H_Sym[[k]]` (would solve list-mode rescale
  perfectly but is a larger model-shape change; separate proposal).
- Adding `T_regularization = "low_rank"` (separate prompt:
  `notes/machima2_t_regularization_low_rank_prompt.md`).
- Changing the `l2` path (it's unaffected by this bug).

## Useful pointers

- AIRE-experiments crash log:
  `logs/machima/run_machima2_brain_joint_Tdense_100000.log`
- AIRE-experiments backup of pre-regularization fit:
  `output/machima2_brain_joint_Tdense_100000.rds.preregularized`
- Original prompt that introduced `frobenius_unit`:
  `notes/machima2_dense_t_regularization_prompt.md` (Option A)
- Diagnostic showing the H_Sym collapse this regularization is meant
  to fix: `output/diagnostic_brain_100000_hsym_summary.csv`
