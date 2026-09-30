# To implement: the neugammazero pdf

Removed from the active pipeline on 2026-09-30. Nothing in the manuscript used
it, so there is no published result to validate against, and `setup/patches/`
now builds an environment without it. Everything needed to put it back is here.

`neugamma` and `gammaneud1` are unaffected and still patched in.

## What it is, and how it differs from neugamma

Both mix a neutral spike below `gammaneu = 1e-4` with a gamma over the rest.
They differ in how the two pieces are combined:

| | neugamma (in use) | neugammazero (removed) |
|---|---|---|
| below 1e-4 | gamma pdf **plus** `pneu/1e-4` | **replaced by** `pneu/1e-4` |
| above 1e-4 | `(1-pneu) x gamma` | `(1-pneu)/pgamma x gamma` |
| `pgamma` | — | `1 - gamma.cdf(1e-4)`, so the gamma part renormalises to 1 over its own support |

So neugamma adds the spike on top of the gamma's own mass below 1e-4, while
neugammazero makes the two pieces disjoint. For a small shape the gamma puts a
lot of mass below 1e-4, so the two are not close.

## Parameters and the ranges used

Three parameters, `shape scale pneu`, the same as neugamma. From
`piperef_stddfe_7_fitgamneug.sh:71-78`; these were the `dadi-cli InferDFE`
`--lbounds`/`--ubounds` and the range the starting-point grid spanned:

```
lower  1e-10 1e-10 0
upper  3     1e8   1
bins   10
```

The only difference from neugamma is the shape ceiling: **3, where gamma and
neugamma get 2**. `bins` is values per parameter in the starting-point grid
(`shape` and `scale` log-spaced, `pneu` linear), before the
`2e-6*Na <= shape*scale <= 2*0.5*Na` plausibility filter.

## Putting it back

1. Add to `config/pipeline.conf` section 3:

```bash
NEUGAMMAZERO_LOWER="${NEUGAMMAZERO_LOWER:-1e-10 1e-10 0}"  # shape scale pneu
NEUGAMMAZERO_UPPER="${NEUGAMMAZERO_UPPER:-3 1e8 1}"
NEUGAMMAZERO_BINS="${NEUGAMMAZERO_BINS:-10}"
```

2. Add `neugammazero` to `PDFS`, and the branch back to `scripts/dfe/run_dfe.sh`:

```bash
  neugammazero)  lower="$NEUGAMMAZERO_LOWER";  upper="$NEUGAMMAZERO_UPPER";  nbins="$NEUGAMMAZERO_BINS" ;;
```

3. Re-add the patch hunks below, and restore the `neugammazero` checks in
   `scripts/checks/test_env.py` (`check_pdfs_registered`,
   `check_neugammazero_is_disjoint`).

## The patch hunks, verbatim

### `dadi/DFE/PDFs.py` — the function

```python
def neugammazero(xx, params):
    """Return Neutral + Gamma distribution"""
    """params = (shape, scale, pneu)
    This is a mixed distribution with a uniform distribution for gamma < 1e-4 with pneu proportion.
    A gamma distribution is used for gamma > 1e-4, with the remaining proportion.
    pneu is the proportion of neutral mutations
    """
    alpha, beta, pneu= params
    xx = np.atleast_1d(xx)
    # Assume gamma < 1e-4 is essentially neutral
    # if gammaneu is not provided, assume it's 1e-4
    # Change it to 2Nas approximately 1e-5, Na is about 18000. 
    gammaneu = 1e-4
    #gammaneu = 0.36
    #gammaneu = 1e-4
    pgamma=1-ssd.gamma.cdf(x=gammaneu, loc=0, a=alpha, scale=beta)
    out = ((1-pneu)/pgamma)*ssd.gamma.pdf(xx,loc=0, a=alpha, scale=beta)
    out[np.logical_and(0 <= xx, xx < gammaneu)] = pneu/gammaneu
    # Reduce xx back to scalar if it's possible
    return np.squeeze(out)
```

Plus the header line it carried: `#2024/11/18 ChenluDi: Added neugammazero
distribution (L41-L51)`.

### `dadi_cli/Pdfs.py` — the four registrations

In `get_dadi_pdf`:

```python
    elif pdf == "neugammazero":
        func = DFE.PDFs.neugammazero
```

In `get_dadi_pdf_params`:

```python
    elif pdf == "neugammazero":
        params = ["shape", "scale", "pneu"]
```

In the name listing:

```python
    print("- neugammazero")
```

In the help text, the string and its branch:

```python
    neugammazero = """
        neugammazero is a mixture of probability density function.
        pneu proportion with selection coefficient <10^-4 and assume to be neutral.
        gamma distribution only for gamma > 1e-4.
        params = [alpha, beta, pneu] = [shape, scale, pneu]
    """
```

```python
    elif pdf_name == "neugammazero":
        print("- neugammazero:\n" + neugammazero)
```
