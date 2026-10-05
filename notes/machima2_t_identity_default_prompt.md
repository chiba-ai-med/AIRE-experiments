# Prompt: Make `fixT=TRUE` (identity T) the default in `Machima2()`

## Audience

This prompt is for another Claude Code session that will modify the
**Machima R package** (https://github.com/kokitsuyuzaki/Machima),
v1.1.0 / SHA `c0612b4`. The user is the package maintainer.

This prompt is intended to be applied **together with**
`notes/machima2_hsym_diagonal_option_prompt.md` (`H_Sym_structure`
option) in a single v1.2.0 release. The two changes together establish
the recommended canonical configuration for cell-type deconvolution of
bulk Hi-C: `fixT = TRUE` + `H_Sym_structure = "diagonal"`.

## Context

`Machima2()` (the symmetric epigenome variant) decomposes paired
scATAC + Hi-C as

```
X_scATAC[k] ≈ W[k] · H_RNA
X_HiC[k]    ≈ G[k] · H_Sym · G[k]^T,    G[k] = T[k] · W[k]
```

Currently `T[k]` is treated as a learned `l_k × n_k` translator between
"ATAC bin space" and "Hi-C bin space". This design is inherited from
the original asymmetric `Machima()`, where the two modalities were
genuinely different feature spaces (gene expression × peak / chromatin
features). For the asymmetric case `T` had a clean interpretation as
"gene → peak translation".

For the **symmetric** Hi-C variant, however, the typical use case
(including AIRE-experiments) re-bins scATAC peaks to the **same
genomic 100 kb grid** as the Hi-C contact matrix, so

```
n_k (= ATAC bin count for chrom k) = l_k (= Hi-C bin count for chrom k)
```

and bin index `i` refers to the **same genomic locus** in both
modalities. There is no "translation" to do — they share the same
coordinate axis. This makes the integration **vertical**
(same cells, shared feature space) rather than diagonal (different
cells, different features). Allowing `T` to be a free non-negative
matrix introduces a learned coordinate transformation between two
spaces that should be identical.

## Symptom (empirical motivation)

In AIRE-experiments / brain 100 kb (Machima v1.1.0, 13 stages × 2 T
variants, completed 2026-05-01 02:07):

### Self-fit reconstruction quality (mean per-chr `rel_frob`, lower better)

| stage              | Identity T | Dense T (frobenius_unit) | Identity better? |
|--------------------|-----------:|-------------------------:|------------------|
| joint              | **0.123**  | **1.000**                | YES (8× better)  |
| jointLowRank       | —          | 0.163                    | (low_rank close) |
| transferFlog       | 0.806      | **81.65**                | YES (catastrophic) |
| transferKraw       | 0.932      | 2.135                    | YES              |
| supervisedGNMF     | 0.933      | 1.047                    | tied             |
| supervisedHfix     | 0.969      | 0.847                    | Dense slightly  |
| supervisedWinit    | 0.938      | 0.738                    | Dense slightly  |

`rel_frob ≥ 1` means the reconstruction is worse than predicting all
zeros. **5 of 7 stages have Dense `rel_frob` ≥ 1 with frobenius_unit**.

The cause: `frobenius_unit` rescales `T` so that
`||T[k]||_F ≈ 1` per chrom, vs Identity's natural `||I_{n_k}||_F = √n_k ≈ 44`
per chrom (for 100 kb on chr1). The constrained Dense T is **44× smaller
in magnitude** than Identity, and `H_Sym` cannot fully compensate
under non-negativity + multiplicative-update constraints. Without
`frobenius_unit` the Dense T fits ~equally to Identity but produces
the H_Sym collapse symptom that motivated v1.1.0 in the first place
(see `notes/machima2_dense_t_regularization_prompt.md`).

### Cell-type discrimination (held-out Bonev sorted-bulk)

`diff_pearson = npc_pearson − cn_pearson`, the per-chrom Pearson
gap between NPC-masked and CN-masked reconstructions vs sorted Bonev
NPC-only / CN-only bulk Hi-C
(`output/sorted_validation_brain_*_100000.csv`):

| stage × T            | diff_pearson | note |
|----------------------|-------------:|------|
| **joint × Tidentity**| **0.951**    | only stage with clear discrimination |
| all 7 Tdense stages (joint, jointLowRank, transfer*, supervised*) | ≈ 0.094 | indistinguishable; all converge to fitting bulk Hi-C without cell-type structure |
| 5 non-joint Tidentity stages | ≈ −0.07 | total failure (negative Pearson) |

The 7 Tdense stages produce **identical sorted-bulk fingerprints**
(`npc_pearson ≈ 0.936, cn_pearson ≈ 0.842`) regardless of Stage A
init strategy. This means the learned T is absorbing all the
modeling freedom and converging to the same bulk-Hi-C fit, irrespective
of how the W factor was initialized. **The stage axis becomes
informationless once Dense T is allowed.**

### Conclusion

Learned T is **strictly harmful** for the symmetric Hi-C case in
AIRE-experiments:
- Hurts reconstruction (Identity fits 8× better in joint)
- Destroys cell-type discrimination (Dense flattens all stages to a
  single bulk fit)
- Justified theoretically by the shared-coordinate structure of
  multiome scATAC + Hi-C

## Goal

Change `Machima2()`'s **default** to `fixT = TRUE` (identity T,
no learning), with auto-construction of identity T when the user
omits the `T` argument. Make learned T explicit opt-in for users
who really need it (e.g., scATAC and Hi-C at different resolutions).

## API surface change

```r
Machima2(
  X_RNA, X_Epi, label = NULL,
  T    = NULL,
  fixT = TRUE,                    # CHANGED: was FALSE in v1.1.0
  ...
)
```

### Behaviour

| `T` arg | `fixT` | new behaviour | back-compat |
|---------|--------|---------------|-------------|
| `NULL` (default) | `TRUE` (new default) | auto-construct `T[k] = diag(l_k)` per chrom (or `diag(l)` in matrix mode), freeze | NEW canonical default |
| `NULL` | `FALSE` (explicit) | random non-negative dense T, learn (= old default) | unchanged behaviour, but user must opt in |
| user-supplied list/matrix | `TRUE` | use as init, freeze | unchanged |
| user-supplied list/matrix | `FALSE` | use as init, learn | unchanged |

### Deprecation warning

When the user opts into learned T (`fixT = FALSE`) for the symmetric
variant **without specifying `T_regularization`**, emit a `warning()`:

```
Learning T (fixT = FALSE) without T_regularization is not recommended
for paired scATAC + Hi-C with shared bin grids. The unconstrained dense
T tends to absorb modeling freedom and degrade cell-type discrimination
(see Machima vignette / NEWS for v1.2.0). Consider:
  - fixT = TRUE (default in v1.2.0+) for shared-grid data
  - T_regularization = "frobenius_unit" or "low_rank" if T must be learned
```

Do **not** error — some users have legitimate use cases (e.g.,
ATAC at 10 kb + Hi-C at 100 kb, where T is a fixed downsampler).

## Implementation

### `R/Machima2.R`

```r
Machima2 <- function(..., fixT = TRUE, ...) {   # changed default
  ...
}
```

### `R/initMachima2.R`

In both `.initMachima2_Matrix` and `.initMachima2_List`, when
`fixT == TRUE && is.null(T)`, construct identity T(s) automatically:

```r
# Inside .initMachima2_List:
if (fixT && is.null(T)) {
  T <- lapply(X_Epi, function(x) diag(nrow(x)))
}
# Inside .initMachima2_Matrix:
if (fixT && is.null(T)) {
  T <- diag(nrow(X_Epi))
}
```

This is the same construction `run_machima2.R` (the AIRE-experiments
caller) currently does manually before calling `Machima2()` — moving
it into Machima itself simplifies callers.

### `R/checkMachima2.R`

Add the deprecation warning:

```r
if (!fixT && (T_regularization == "none") && is.null(T)) {
  warning(
    "Learning T (fixT = FALSE) without T_regularization is not ",
    "recommended for paired scATAC + Hi-C with shared bin grids. ",
    "The unconstrained dense T tends to absorb modeling freedom ",
    "and degrade cell-type discrimination. Consider fixT = TRUE ",
    "(default in v1.2.0+) or T_regularization = ",
    "\"frobenius_unit\" / \"low_rank\"."
  )
}
```

### Roxygen on `Machima2()`

Update `@param fixT`:

```
@param fixT Logical. If TRUE (default in v1.2.0+), T is fixed during
iteration: when T is also NULL, an identity matrix is auto-constructed
per chrom. If FALSE, T is learned as a non-negative dense matrix.
For paired scATAC + Hi-C on a shared bin grid, the recommended setting
is fixT = TRUE; the learned-T mode is intended for resolution-mismatched
or asymmetric cases. (Default: TRUE)
```

Update `@details`:

```
The symmetric variant of Machima2 assumes scATAC features and Hi-C bins
share the same genomic coordinate system (typical for paired multiome
where ATAC peaks are re-binned to the Hi-C bin grid). Under this
assumption, T is naturally the identity matrix, and learning T as a
free non-negative matrix introduces a fictitious coordinate
transformation that degrades both reconstruction quality and cell-type
discrimination. The v1.2.0 default is therefore fixT = TRUE.

For asymmetric or resolution-mismatched cases (e.g., scATAC at 10 kb
binned vs Hi-C at 100 kb binned), fixT = FALSE with a non-trivial
T_regularization remains supported.
```

## NEWS entry (v1.2.0)

```
VERSION 1.2.0
------------------------
   o BREAKING: Machima2() default for fixT changed from FALSE to TRUE.
     Auto-constructs identity T per chrom when T = NULL && fixT = TRUE.
     Existing callers that rely on learned T must explicitly set
     fixT = FALSE.

   o Deprecation warning added when fixT = FALSE is used without
     T_regularization for paired-modality cases.

   o (Companion change) Added H_Sym_structure = c("symmetric", "diagonal")
     option for opt-in symmetric CP-style decomposition.
```

## Acceptance test

Add `tests/testthat/test-Machima2-default-fixT.R`:

1. **New default**: `Machima2(X_RNA, X_Epi)` runs without explicit
   `fixT` and produces `res$T` equal to identity (matrix mode) or
   list of identities (list mode). `res$.meta$fixT == TRUE` (or
   equivalent inspection).
2. **Auto-identity construction**: `Machima2(X_RNA, X_Epi, T = NULL,
   fixT = TRUE)` returns `res$T` as identity without the caller
   passing it.
3. **Backward compat (explicit learned T)**: `Machima2(X_RNA, X_Epi,
   fixT = FALSE)` runs with the old behaviour (random dense T,
   learned), produces identical results to v1.1.0 modulo random
   seed (within numerical tolerance).
4. **Deprecation warning**: `expect_warning(Machima2(X_RNA, X_Epi,
   fixT = FALSE), "fixT = FALSE.*not recommended")`.
5. **No warning when regularization is set**: `Machima2(X_RNA, X_Epi,
   fixT = FALSE, T_regularization = "low_rank", T_rank = 4)` runs
   silently (no deprecation warning), since the user has explicitly
   selected a regularization strategy.
6. **AIRE-experiments brain 100 kb regression**: with new defaults,
   the joint stage reproduces the existing `joint × Tidentity`
   numbers (rel_frob ≈ 0.12, diff_pearson ≈ 0.95).

## Migration guidance

Downstream callers fall into three categories:

1. **Callers that already set `fixT = TRUE` explicitly** (e.g., the
   AIRE-experiments `run_machima2.R` `T_variant = "identity"` branch):
   no change needed; will continue to work.
2. **Callers that omit `fixT`** and were getting learned T by default:
   will silently switch to identity T. This is the intended
   behavioural change. Documented in NEWS.
3. **Callers that set `fixT = FALSE` explicitly**: continue to get
   learned T, plus a deprecation warning if no `T_regularization`
   is set.

`run_machima2.R` in AIRE-experiments will need a minor cleanup to
drop the manual identity construction (Machima now does it), but the
semantic behaviour is unchanged.

## Out of scope

- Removing `T`, `fixT`, `T_regularization`, or `T_rank` from the API.
  All remain available for advanced users and resolution-mismatched
  cases.
- Changing `Machima()` (asymmetric variant) defaults. The asymmetric
  case has genuinely different feature spaces (e.g., gene × peak),
  so learned T is biologically meaningful there.
- Changing `H_Sym` defaults — that's the companion prompt
  (`notes/machima2_hsym_diagonal_option_prompt.md`).
- Modifying the `T_regularization = "low_rank"` or
  `T_regularization = "frobenius_unit"` paths themselves.

## Useful pointers

- AIRE-experiments self-fit metrics:
  `output/eval_brain_*/metrics.csv` (column `mean_rel_frob_per_chr`)
- AIRE-experiments held-out validation:
  `output/sorted_validation_brain_*_100000.csv`
- AIRE-experiments H_Sym diag/off-diag analysis:
  `output/diagnostic_brain_100000_hsym_diag_ratio.csv`
- Companion spec change (apply together):
  `notes/machima2_hsym_diagonal_option_prompt.md`
- Math reference for the symmetric NMTF model:
  `notes/machima2_update_rules.tex` / `.pdf`
- Earlier T-related prompts (already landed in v1.0.0/v1.1.0):
  - `notes/machima2_dense_t_regularization_prompt.md`
  - `notes/machima2_frobenius_unit_listmode_fix_prompt.md`
  - `notes/machima2_t_regularization_low_rank_prompt.md`

These earlier prompts solved the H_Sym collapse symptom of learned T
mathematically (rescale / regularize), but the present empirical
results show that for paired same-grid data the cleanest fix is to
not learn T at all. v1.2.0 makes that the default while keeping the
earlier machinery available for cases that need it.
