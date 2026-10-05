# Prompt: Add T-regularization to Machima2 (Dense T variant)

## Audience

This prompt is for another Claude Code session that will modify the
**Machima2 R package** (https://github.com/kokitsuyuzaki/Machima).
The user is the package maintainer.

## Context

Machima2 model:

```
X_RNA[k]  ~=  W[k] . H_RNA
X_Epi[k]  ~=  (T[k] . W[k]) . H_Sym . (T[k] . W[k])^T
```

`T[k]` is `l_k x n_k`. Two variants are currently supported:

- `T_variant = "identity"` — `T[k] = diag(l_k)` (frozen, `fixT = TRUE`)
- `T_variant = "dense"`    — `T[k]` is a free dense matrix, learned

Used in: AIRE-experiments / brain 100 kb shake-down, both variants run for
6 NMF stages each (12 total).

## Symptom (problem to solve)

For **Dense T** variants, after fitting, `H_Sym` collapses to near-zero
absolute values (~ 1e-12), while `T[k]` becomes very large. The product
`(T.W) . H_Sym . (T.W)^T` is well-scaled (Hi-C reconstruction `rel_frob`
~= 0.12 — good fit), but **all the structure is in T, not in H_Sym**.

Concrete evidence (from `output/diagnostic_brain_100000_hsym_summary.csv`
in the AIRE-experiments repo):

| stage             | T        | diag_mean (npc) | diag_mean (cn) | cross-within   |
|-------------------|----------|-----------------|----------------|----------------|
| supervisedHfix    | dense    | 9.83e-13        | 9.48e-13       | 3.36e-13       |
| supervisedHfix    | identity | 7.05e-08        | 7.59e-15       | 2.21e-11       |
| joint             | dense    | 1.20e-03        | 1.16e-03       | -3.06e-04      |
| joint             | identity | 479             | 41.5           | 34.9           |

Identity T runs have non-trivial H_Sym structure (several orders of
magnitude). Dense T runs have H_Sym essentially at numerical noise scale.

Downstream effect: differential Pearson against held-out per-celltype
Bonev bulks is ~= 0 across all 12 stages. Cell-type-specific contact
patterns are not encoded in `H_Sym`, so component-mask reconstructions
`(T.W) . M_g . H_Sym . M_g^T . (T.W)^T` carry no cell-type discrimination.

## Why this happens

`T` has `l_k * n_k` (~ millions per chrom at 100 kb) free parameters.
`H_Sym` has `J^2` (~ 49) free parameters. The optimizer routes all
structure to T because (a) it has vastly more capacity, (b) the loss
landscape has a flat direction `T <- T*c, H_Sym <- H_Sym/c^2` for any
`c > 0`, and the optimizer drifts toward `c = large` over iterations.

## Required change

Add an option to **constrain T's scale during fitting** so H_Sym retains
informative magnitude. Implement at least one of these (in order of
preference); allow user to opt in.

### Option A (preferred): per-iteration Frobenius normalization

After each T update step, rescale T per chromosome:

```
T[k] <- T[k] / ||T[k]||_F
H_Sym <- H_Sym * (mean_k ||T[k]||_F)^2   # absorb magnitude into H_Sym
```

(or any equivalent rescale that preserves `T[k] . W[k] . H_Sym . W[k]^T . T[k]^T`
exactly; the point is to fix T's scale and let H_Sym carry the magnitude.)

Pros: zero hyperparameter, no extra loss term, model output unchanged.
Cons: can momentarily destabilize convergence if applied too aggressively.

### Option B: L2 penalty on T

Add `lambda_T * sum_k ||T[k]||_F^2` to the loss. Update rule for T includes
the derivative `2 * lambda_T * T[k]`. User-supplied `lambda_T`; default
should be a small fraction of `||X_Epi||_F^2 / (sum_k l_k * n_k)`.

Pros: standard regularization. Cons: needs hyperparameter tuning.

### Option C (most invasive): low-rank T

Parametrize `T[k] = U[k] V[k]^T` with `U[k]: l_k x r`, `V[k]: r x n_k`,
`r << min(l_k, n_k)`. User-supplied `T_rank`.

Pros: bounded capacity; principled. Cons: rewrites T update.

## Suggested API surface

Add to `Machima2()`:

```r
Machima2(
  ...,
  T_regularization = c("none", "frobenius_unit", "l2", "low_rank"),
  lambda_T         = 0,            # used when T_regularization = "l2"
  T_rank           = NULL,         # used when T_regularization = "low_rank"
  ...
)
```

Defaults must be backward compatible: `T_regularization = "none"`.
Identity T (`fixT = TRUE`) ignores `T_regularization`.

## Acceptance test

The user will re-run AIRE-experiments / brain 100 kb with
`T_regularization = "frobenius_unit"` for the 6 Dense T stages and check:

- `H_Sym` diagonal values are O(1) (not O(1e-12))
- `H_Sym` shows visible NPC-vs-CN structure: `diag_mean(npc) - diag_mean(cn)`
  is non-trivial relative to the overall scale
- Differential Pearson against held-out Bonev sorted bulks moves
  meaningfully off zero (target: > 0.2 for at least one stage)
- Hi-C `rel_frob` does not regress significantly (stays around 0.12-0.20)

## Out of scope

- Changing the model's `X_RNA = W H_RNA` side
- Changing the `init_*` API (already added in a previous round)
- Changing fixT / Identity T path

## Useful pointers

- Existing `init_*` API was added per `notes/machima2_init_api_prompt.md`
- AIRE-experiments diagnostic that produced the H_Sym evidence:
  `src/diagnose_hsym_and_bonev_diff.R`
- The shake-down PDF (Japanese) summarising Dense T behavior:
  `output/summary_brain_100000.pdf` page 9 ("考察と次の検討事項")
