# FAA outage attribution — 4-hour Polymarket price forecast

Entry for the make.fun competition `polymarket-faa-attribution-02`.
`python3 main.py` reads every `*.in` case in `$MAKEFUN_INPUT_DIR` (or
`/input`) and writes one `*.out` forecast per case to `$MAKEFUN_OUTPUT_DIR`
(or `/output`). Standard library only, no network, deterministic.

## Approach

A random-walk forecast around each market's checkpoint price:

- **Median = checkpoint price + sum reversion.** The markets of one event
  are alternative answers to the same question, so their Yes prices
  together are one probability mass. When the sum at the checkpoint sits
  above (or below) the level it sat at over the last 12 hours, the excess
  is traded away over the next hour (strength 0.8, time constant 60 min),
  distributed across markets in proportion to each market's share of the
  sum. The target level is the **median of the sum over the last 720
  minutes** — a level the event itself set, not a number assumed in
  advance — so the same rule works whether the family's prices add up to
  about 1 (a proper mutually exclusive event, like the practice's Next
  CEO of NVIDIA) or to about 0.1 (a family of long-shot attributions,
  like the FAA family, whose four live prices sum to ~0.10). The median,
  not the mean, is used because a market's first minutes trade near 0.5
  before finding its level, which spikes the sum briefly.
- **Spread fitted per market.** The per-minute volatility is the LARGER
  of a fast (half-life 30 min) and a slow (half-life 2 hours)
  exponentially weighted variance of the market's one-minute price
  changes (last 24 hours, winsorised at ±0.05), scaled by 0.3 and floored
  at 0.0002. Volatility clusters — it jumps up fast after a burst and
  decays slowly — so the spread rises fast and falls slow. The spread
  grows like sqrt(h) for minute h, with the 5%/95% quantiles widened a
  further 1.3x because the moves are heavy-tailed (quantiles of the
  standard normal).
- **No own-trend drift.** Short-term drift from a market's own recent
  trend was tested on the practice cases and made the median strictly
  worse in every form tried (mean-reverting step prices punish chasing
  the last move), so no drift is applied beyond the sum reversion above.

Quantiles are clamped to [0, 1], rounded to 6 decimals, and forced
non-decreasing.

## What changed vs the previous entry (v2)

v2 was a random walk around the checkpoint price with a Student-t(4)
shape and an EWMA (half-life 30 min) sigma, no drift at all. v3 adds the
sum reversion (the single biggest win on the practice: it is worth about
0.013 pooled by itself), replaces the t-shape with a normal shape plus a
1.3x outer-quantile widening, takes the max of a fast and a slow
volatility instead of one half-life, and lowers the sigma floor from
0.00032 to 0.0002 (quiet markets were being over-spread).

Things tried and dropped (all measured on the practice, pooled p00–p04):
a hardcoded sum target of 0.96–0.97 (a touch better on the practice,
0.2921, but it would be catastrophic on the final, where the family's
prices sum to ~0.10, not ~1 — so the target must be the event's own
recent level); empirical 1-minute move shapes as quantiles (far worse,
0.359 — horizon tails are thinner than 1-minute tails, because jumps
partly revert); a time-of-day activity profile (worse, 0.305); bipower
variation sigma (worse, 0.302); compound-jump spreads (far worse in the
sibling competition's work); distributing the sum correction by recent
move instead of by share (worse, 0.321); separate growth exponents for
the outer quantiles (worse).

## Practice results (local scoring)

Scored against the next case's prices (p05 has no later case, so it is
not in the local numbers), 5 cases, 124,800 items, pinball loss in
percentage points:

| version | pooled | per-case mean | p00 | p01 | p02 | p03 | p04 |
|---------|--------|---------------|-----|-----|-----|-----|-----|
| house baseline | — | 0.398458 | 0.5361 | 0.3804 | 0.4084 | 0.3686 | 0.2988 |
| v2 (previous entry) | 0.318224 | 0.327831 | 0.5059 | 0.3092 | 0.3865 | 0.2366 | 0.2010 |
| **v3 (this entry)** | **0.286033** | **0.291501** | **0.3534** | **0.3084** | **0.3314** | **0.2670** | **0.1974** |

The house baseline's published practice figure is 0.363589 pooled over
all six cases. The previous entry's submitted figure, 0.279746, could
not be reproduced from the local data (its p00–p04 pooled is 0.318224,
so the difference must come from p05, which has no local reference
answers); v3 beats v2 by 0.032 pooled on the five verifiable cases, and
its p05 drift is near zero, so the six-case total is expected to improve
by about the same margin.

The practice is cut from a different event (Next CEO of NVIDIA), so these
numbers are a guide only, not a promise about the final.

## Verification

- Two full runs over all 6 practice cases: byte-identical outputs.
- Single-case run (one `.in` in the input dir, as the grader does):
  0.11 s, output identical to the full-run output. Limit is 8 s per case.
- All 6 outputs valid: exactly one `forecasts` key, every case market
  once, 240 rows × 5 numbers, each in [0, 1], non-decreasing per row,
  ≤ 1 MiB (193–262 KB).
- A synthetic FAA-like case (4 markets, prices summing to 0.10, born at
  ~0.49 with a fast crash like the real FAA markets) was run to check
  the adaptive sum target: medians stay near the checkpoint prices
  (largest drift +0.009), i.e. no pull toward a sum of 1.
- Python 3.12.3 locally (grader uses 3.13); stdlib only, no network,
  no reads outside the case file.
