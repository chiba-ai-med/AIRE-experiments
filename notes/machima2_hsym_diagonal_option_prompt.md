# Prompt: Add `H_Sym_structure` option to constrain `H_Sym` to diagonal

## Audience

This prompt is for another Claude Code session that will modify the
**Machima R package** (https://github.com/kokitsuyuzaki/Machima),
v1.1.0 / SHA `c0612b4`. The user is the package maintainer.

## Context

`Machima2()` (the symmetric epigenome variant) decomposes a Hi-C contact
matrix as

```
X_HiC[k] ≈ G[k] · H_Sym · G[k]^T,    G[k] = T[k] · W[k]
```

where `H_Sym` is `J × J` symmetric non-negative. Currently `H_Sym` has
`J(J+1)/2` free parameters; off-diagonal entries `H_Sym[i, j]` (`i ≠ j`)
parameterize cross-cell-type contacts.

This is mathematically a **symmetric Tucker** core: each off-diagonal
entry contributes a rank-1 non-symmetric piece `g_i g_j^T`, paired with
its transpose `g_j g_i^T` to form a symmetric rank-2 block representing
"cell-type i ↔ cell-type j cross-contact propensity".

## Symptom (motivation for the change)

In AIRE-experiments / brain 100 kb (after Machima v1.1.0 frobenius_unit
fix + 7 stages × 2 T variants + 1 jointLowRank Tdense = 13 fits), the
ratio `||off-diag(H_Sym)||_F / ||diag(H_Sym)||_F` was computed across
all stages (`output/diagnostic_brain_100000_hsym_diag_ratio.csv`):

- **10 / 13** stages have ratio ≥ 0.5 (off-diagonal substantial or
  dominant)
- **3 / 13** Tdense supervised stages have ratio < 0.1 — but inspection
  shows these are **single-component collapses** (one diagonal entry
  carries ~99% of the trace), not healthy diagonal-dominance

**Biological interpretation problem.** For cell-type deconvolution of
bulk Hi-C, each leiden cluster is a putative cell type. The "interaction"
between cell-type `i` and cell-type `j` (off-diagonal entry of `H_Sym`)
is hard to assign a clean biological meaning:

- It is not a contact "between two cells of different types" (single-cell
  Hi-C contacts are within one cell)
- It would have to represent a **co-occurrence pattern** of contacts
  attributed to type `i` at locus `a` and type `j` at locus `b` — but
  the bulk Hi-C input is a sum over cells, so this "co-occurrence" is
  ambiguous between (a) physically interacting cell types in tissue
  architecture and (b) a numerical artifact of overlapping component
  loadings

The user has decided that for this application (cell-type deconvolution
of bulk Hi-C), the **per-component independent contribution** model
(diagonal `H_Sym`) is biologically more defensible than the cross-type
core (full symmetric `H_Sym`). This corresponds to constraining the
factorization to **symmetric CP** (a.k.a. INDSCAL, sym-PARAFAC) form:

```
X_HiC[k] ≈ Σ_i h_i · g_i[k] · g_i[k]^T
```

instead of the current symmetric Tucker form

```
X_HiC[k] ≈ Σ_{i,j} H_Sym[i,j] · g_i[k] · g_j[k]^T.
```

## Goal

Add a new argument **`H_Sym_structure`** to `Machima2()` that lets the
user constrain `H_Sym` to be diagonal-only. Default keeps current
behaviour (back-compat).

## API surface

```r
Machima2(
  ...,
  H_Sym_structure = c("symmetric", "diagonal"),
  ...
)
```

Behaviour:

- `"symmetric"` (default): existing behaviour. `H_Sym` is `J × J`
  symmetric non-negative. No change.
- `"diagonal"`: `H_Sym` is constrained to `diag(h_1, ..., h_J)`. All
  off-diagonal entries are kept identically zero throughout the
  iteration. This makes the model symmetric CP / sym-PARAFAC.

`H_Sym_structure = "diagonal"` is **independent of and orthogonal to**
all other knobs (`fixT`, `T_regularization`, `init_*`, `Beta`, `J`,
horizontal mode, etc.). Specifically:

- Compatible with `fixT = TRUE` (identity T) — most natural pairing
- Compatible with `T_regularization = "frobenius_unit"` — the rescale
  `H_Sym ← c^2 · H_Sym` preserves diagonal structure
- Compatible with `T_regularization = "low_rank"` — diagonal `H_Sym`
  + low-rank T is a valid sym-CP-with-low-rank-G factorization
- Compatible with `init_H_Sym`: if user supplies non-diagonal
  `init_H_Sym`, the off-diagonal is silently projected to zero with a
  warning (preferred over an error — supports the workflow of warm-
  starting from a previous symmetric fit)

## Implementation

### Recommended approach: project after the standard MU step

In `R/Machima2.R` main loop, after the existing `H_Sym` MU update +
symmetrize step:

```r
if(!fixH_Sym){
    H_Sym <- .updateH_Sym(...)                    # existing
    H_Sym <- (H_Sym + t(H_Sym)) / 2               # existing
    if(H_Sym_structure == "diagonal"){
        H_Sym <- diag(diag(H_Sym))                # NEW
    }
}
```

This is the simplest, lowest-risk implementation: 3 added lines, reuses
the existing update code, the projection step is idempotent and fast
(`O(J^2)` per iter, negligible).

### Alternative (more efficient but more code): diagonal-only update

When `H_Sym_structure == "diagonal"`, the MU update simplifies because
only the J diagonal entries need to be updated. The numerator and
denominator collapse to scalar quantities per diagonal index `i`:

```
numer_i = sum_k pi_HiC[k] * (g_i[k]^T (X_HiC[k] ⊙ Y_HiC[k]^(β-2)) g_i[k])
denom_i = sum_k pi_HiC[k] * (g_i[k]^T Y_HiC[k]^(β-1) g_i[k]) + L_terms
h_i ← h_i * (numer_i / denom_i)^ρ
```

This avoids forming the full `J × J` numerator/denominator matrices and
is asymptotically `O(J)` updates instead of `O(J^2)`. For `J = 7` this
saves nothing measurable; for large `J` it would matter.

**Choose the project-after approach** for v1.2.0 (less code, easier
to review, easier to test). Optimize to diagonal-only updates only if
profiling shows it's a bottleneck.

### Initialization

In `R/initMachima2.R`:

- The existing init code constructs `H_Sym` (random or from
  `init_H_Sym` user arg). At the end of the init function, add:

```r
if(H_Sym_structure == "diagonal"){
    if(!is.null(init_H_Sym) && any(abs(init_H_Sym - diag(diag(init_H_Sym))) > 1e-12)){
        warning("init_H_Sym has non-zero off-diagonal entries; ",
                "projecting to diagonal because H_Sym_structure='diagonal'")
    }
    H_Sym <- diag(diag(H_Sym))
}
```

### Validation

In `R/checkMachima2.R`:

```r
H_Sym_structure <- match.arg(H_Sym_structure)   # in Machima2()
# In .checkMachima2():
stopifnot(H_Sym_structure %in% c("symmetric", "diagonal"))
```

## Documentation

### Roxygen on `Machima2()`

Add to `@param`:

```
@param H_Sym_structure Constrain the J×J symmetric core H_Sym to be
"symmetric" (default; full J(J+1)/2 free parameters; corresponds to
symmetric Tucker) or "diagonal" (only J diagonal entries; corresponds
to symmetric CP / sym-PARAFAC). The diagonal option is recommended for
cell-type deconvolution of bulk Hi-C, where off-diagonal entries are
hard to assign a clean biological meaning. Compatible with all other
flags. (Default: "symmetric" for backward compatibility.)
```

Add to `@details`:

```
The diagonal mode (H_Sym_structure = "diagonal") implements a
symmetric CP-like factorization

  X_Epi[k] ≈ sum_i h_i * g_i[k] * g_i[k]^T

where each component i has an independent contact pattern g_i[k]
weighted by h_i = H_Sym[i,i]. This is the natural choice when
components are interpreted as independent cell types and cross-type
"interaction" is not biologically meaningful for the application.

The default symmetric mode allows off-diagonal H_Sym[i,j] to capture
co-occurrence of contacts between components i and j, useful when
the dataset has spatially co-localized cell types whose contacts
must be modeled jointly.
```

### NEWS entry

```
VERSION 1.2.0
------------------------
   o Added H_Sym_structure = c("symmetric", "diagonal") to Machima2()
     to allow constraining H_Sym to a diagonal core (sym-CP /
     sym-PARAFAC), which is more biologically defensible for
     cell-type deconvolution of bulk Hi-C
```

## Acceptance test

Add `tests/testthat/test-Machima2-Hstruct.R`:

1. **Diagonal stays diagonal**: `Machima2(..., H_Sym_structure = "diagonal", num.iter = 30)`
   on a small random fixture; verify `all(res$H_Sym == diag(diag(res$H_Sym)))`
   exactly (not approximately — the projection makes them identically zero).
2. **Backward compat**: `Machima2(..., H_Sym_structure = "symmetric")`
   with a fixed seed produces identical output to the existing
   `Machima2(...)` call without the new arg (modulo random init).
3. **Compatible with fixT**: `Machima2(..., fixT = TRUE,
   H_Sym_structure = "diagonal")` runs and returns a diagonal H_Sym.
4. **Compatible with frobenius_unit**: ditto with
   `T_regularization = "frobenius_unit"`; rescale should not introduce
   off-diagonal mass.
5. **Compatible with low_rank**: ditto with
   `T_regularization = "low_rank", T_rank = 2*J`.
6. **Init projection warning**: passing
   `init_H_Sym = matrix(runif(J*J), J, J)` (asymmetric, full)
   together with `H_Sym_structure = "diagonal"` triggers the warning
   and the run completes with a diagonal `H_Sym`.
7. **Reconstruction sanity**: on the AIRE-experiments brain 100 kb
   joint stage with `H_Sym_structure = "diagonal"`, mean per-chr
   `rel_frob` of the Hi-C reconstruction is within a factor of 2 of
   the symmetric variant (diagonal has less capacity, but should not
   collapse).

## Empirical context for the acceptance test

The AIRE-experiments brain 100 kb / joint × Tidentity benchmark gives
the following held-out Bonev sorted-bulk Pearson with the **current
symmetric `H_Sym`** (from `output/sorted_validation_brain_joint_Tidentity_100000.csv`):

```
                       npc_pearson   cn_pearson   diff_pearson
joint × Tidentity        0.911       -0.040         0.951
```

Re-running this stage with `H_Sym_structure = "diagonal"` should:

- Keep `npc_pearson > 0.7` (don't lose NPC reconstruction)
- Either improve or not regress `cn_pearson` (currently bad at -0.04)
- Hopefully **improve `diff_pearson > 0` while having both NPC and
  CN above zero** — i.e., differentiate cell types via independent
  diagonal weights rather than via off-diagonal cross-talk.

This is a "would have been the right model all along" kind of
acceptance test; it's not strict, but if diagonal is much worse than
symmetric, the change is suspect.

## Out of scope

- Changing the **default** to "diagonal" (this prompt keeps default
  "symmetric" for back-compat; default change is a separate decision
  after empirical validation on more datasets)
- Block-diagonal `H_Sym` or low-rank `H_Sym` (tier-2 future work)
- Refactoring the `H_Sym` update math itself
- The `Machima` (asymmetric, non-`H_Sym` variant) function — it has
  no `H_Sym` argument

## Useful pointers

- AIRE-experiments diagnostic CSV with the diag/off-diag ratio:
  `output/diagnostic_brain_100000_hsym_diag_ratio.csv`
- AIRE-experiments held-out validation CSVs:
  `output/sorted_validation_brain_*_100000.csv`
- Discussion of the sym-Tucker vs sym-CP tradeoff:
  `notes/machima2_update_rules.tex` / `notes/machima2_update_rules.pdf`
  Section "Update for H_Sym" (frobenius_unit interaction is in the
  same section)
- Related existing prompts:
  - `notes/machima2_frobenius_unit_listmode_fix_prompt.md` (compatible)
  - `notes/machima2_t_regularization_low_rank_prompt.md` (compatible)
  - `notes/machima2_hsym_dimnames_na_prompt.md` (orthogonal: dimnames
    is independent of structure)
- The symmetric NMTF / sym-PARAFAC literature: Hoyer 2004 ("Non-negative
  Matrix Factorization with Sparseness Constraints"), Bro 1997
  (PARAFAC review), Comon et al. 2008 (symmetric tensor decomposition)
