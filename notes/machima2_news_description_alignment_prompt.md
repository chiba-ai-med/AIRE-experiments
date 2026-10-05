# Prompt: Align `inst/NEWS` and `DESCRIPTION` version history

## Audience

This prompt is for another Claude Code session that will modify the
**Machima R package** (https://github.com/kokitsuyuzaki/Machima),
v1.0.0 / SHA `c1d1804`. The user is the package maintainer.

## Context

In the recent two-commit run that landed `init_*` (commit `e69cd23`)
and `T_regularization` + NA dimnames (commit `c1d1804`), the
`DESCRIPTION` and `inst/NEWS` files drifted out of sync.

`DESCRIPTION` (HEAD `c1d1804`):

```
Package: Machima
Version: 1.0.0
```

`inst/NEWS` (HEAD `c1d1804`):

```
VERSION 1.0.0
------------------------
   o Added T_regularization param to Machima2()
     ("frobenius_unit", "l2") to prevent H_Sym collapse
     when T is a free dense matrix
   o Fixed NA in dimnames of H_Sym/H_RNA when J > num_celltypes
     (falls back to "comp_j")

VERSION 0.99.1
------------------------
   o Added Machima2() for symmetric epigenomic matrices
   o Added init_W_RNA/init_H_RNA/init_H_Sym, fixH_Sym to Machima2()
   ...
```

So `NEWS` records two distinct releases (`0.99.1` for the init API,
`1.0.0` for T_regularization), but `DESCRIPTION` only shows `1.0.0`.
The `e69cd23` commit (which `NEWS` calls `0.99.1`) actually shipped
with `DESCRIPTION: Version: 1.0.0` already — verify by
`git show e69cd23:DESCRIPTION`. If that's the case, `0.99.1` was never
actually a released version: the bump-to-1.0.0 happened in `e69cd23`
but `NEWS` only documented it in `c1d1804`.

## Goal

Make the version history coherent with the actual git tags / commits.
Two acceptable resolutions — pick whichever matches the user's
intent for downstream `remotes::install_github` consumers:

### Resolution A (preferred): `0.99.1` is real, `1.0.0` is the T_reg release

- Amend `e69cd23`'s `DESCRIPTION` to set `Version: 0.99.1` (so the
  init-API release ships under that version).
- Keep `c1d1804`'s `DESCRIPTION` at `Version: 1.0.0`.
- `NEWS` is already correct under this resolution; no change.
- Tag commits: `git tag v0.99.1 e69cd23` and `git tag v1.0.0 c1d1804`,
  push tags.

### Resolution B: `0.99.1` was never a release, fold into `1.0.0`

- `DESCRIPTION` already says `1.0.0` from `e69cd23`; no change.
- Edit `inst/NEWS` to merge the `VERSION 0.99.1` and `VERSION 1.0.0`
  blocks into a single `VERSION 1.0.0` entry containing all the
  bullets from both.
- Tag only `git tag v1.0.0 c1d1804`.

## Recommended choice

Resolution A. Reasons:

- Bioconductor / CRAN convention: every commit that bumps `DESCRIPTION`
  is a release; every release gets a `NEWS` entry with the matching
  version header. The init API was a substantial public-API addition
  and deserves its own release row in the NEWS for downstream
  archeology.
- `remotes::install_github("kokitsuyuzaki/Machima@v0.99.1")` becomes a
  meaningful pin for users who want the init API but not the T_reg
  changes.
- AIRE-experiments' own `setup_r_packages.R` does
  `remotes::install_github("kokitsuyuzaki/Machima")` (HEAD-tracking),
  so neither resolution affects it operationally — this is purely
  about future-archeology hygiene.

## Required changes (under Resolution A)

1. **Rewrite the `e69cd23` commit's `DESCRIPTION`** so it carries
   `Version: 0.99.1` instead of `1.0.0`. Mechanically:

   ```sh
   git checkout e69cd23 -- DESCRIPTION
   sed -i 's/^Version: 1\.0\.0$/Version: 0.99.1/' DESCRIPTION
   git commit --amend --no-edit DESCRIPTION   # only if e69cd23 is HEAD~1
   ```

   In practice, `e69cd23` is not HEAD, so do this via
   `git rebase -i 6cd73ae` (interactive) and `edit` the e69cd23
   commit. Or, more safely: leave history alone and just push tags
   pointing at the existing commits, and add a single new commit
   correcting only the NEWS (Resolution B-lite — see below).

2. **`git tag v0.99.1 e69cd23 && git tag v1.0.0 c1d1804`** then
   `git push --tags`.

3. **Add a `VERSION 1.0.0`-block release-date line** to `inst/NEWS`
   (currently the headers are versionless on the date side):

   ```
   VERSION 1.0.0    (2026-04-30)
   ------------------------
   ...
   VERSION 0.99.1   (2026-04-29)
   ------------------------
   ...
   ```

## Out of scope

- Reorganising `inst/NEWS` into NEWS.md (current file is plain `NEWS`,
  matches package's existing format)
- Changing the version-bump cadence policy
- Adding a CHANGELOG/CONTRIBUTING

## Acceptance criteria

1. `git tag --list 'v*'` returns at least `v0.99.1` and `v1.0.0`.
2. `git show v0.99.1:DESCRIPTION | grep ^Version` returns `Version: 0.99.1`.
3. `git show v1.0.0:DESCRIPTION  | grep ^Version` returns `Version: 1.0.0`.
4. `inst/NEWS` `VERSION` headers match the version field in
   `DESCRIPTION` at that tag.
5. `R CMD check` clean, no new NOTEs.

## Useful pointers

- Two recent commits: `e69cd23` (init API), `c1d1804` (T_reg + NA dimnames)
- Pre-existing baseline: `6cd73ae Machima2 (Symmetric)`
