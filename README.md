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
  grows like sqrt(240·(1 − exp(−h/240))) for minute h (saturating, not
  sqrt(h): prices mean-revert inside the event, mass swapping between
  markets, so long-horizon spread grows slower than a random walk's),
  with the 95% quantile widened a further 1.8x and the 5% quantile not
  widened at all — the moves are heavy-tailed upward only, the 0
  boundary cutting the lower tail off (quantiles of the standard normal
  otherwise).
- **Leader-decay median shift (new in v4).** A market's future price
  distribution is right-skewed (bounded below by 0, occasional jumps
  up), so its conditional median sits below its checkpoint price by an
  amount scaling with its own volatility. The median is pulled down by
  0.5 · w · sigma · grow(h), where w ramps from 0 at half the
  field-average price (sum / number of markets) to 1 at the field
  average. The weight is relative to the event's own field average, so
  it is scale-free and fires for the larger markets of any family,
  whether its prices sum to ~1 or ~0.1.
- **No own-trend drift.** Short-term drift from a market's own recent
  trend was tested on the practice cases and made the median strictly
  worse in every form tried (mean-reverting step prices punish chasing
  the last move), so no drift is applied beyond the sum reversion and
  the leader-decay shift above.

Quantiles are clamped to [0, 1], rounded to 6 decimals, and forced
non-decreasing.

## What changed vs the previous entry (v3)

v3 was a random walk around the checkpoint price with the sum reversion
and a normal shape with both outer quantiles widened 1.3x, spread
growing like sqrt(h). v4 adds three things:

1. **Leader-decay median shift** (the biggest win, worth about 0.008
   pooled by itself): the median of a larger market is pulled down by
   0.5 · w · sigma · grow(h), w ramping 0→1 between half the field
   average and the field average. Found from the loss decomposition:
   the median is the dominant loss term (oracle-median loss ~0.09
   pooled vs 0.286 actual), and the large markets' errors are
   dominated by big downward moves (crashes) that a symmetric random
   walk cannot anticipate. A per-market oracle scan confirmed the
   optimal median shift is negative for the markets that later crashed
   and near zero for the small field.
2. **Saturating spread growth**: sqrt(240·(1 − exp(−h/240))) instead of
   sqrt(h) (worth about 0.001 pooled).
3. **Asymmetric outer tails**: 95% widened 1.8x, 5% not widened (was
   1.3x both) — the lower tail is cut off by the 0 boundary, the upper
   tail is long (worth about 0.001 pooled).

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
the outer quantiles (worse); **new in v4's search:** median drift toward
the market's own trailing median (mixed: helped p01/p04, hurt p03 badly,
best 0.2854); drift toward the event's target sum times the market's
recent share (far worse, 0.377–0.584 — the recent share is polluted by
new markets' first minutes near 0.5); drift toward an equal share of the
sum (far worse, 0.319–0.493); lognormal (multiplicative) quantile shape
(worse, 0.295–0.376); trailing-median/trailing-mean anchor instead of
the last price (far worse, 0.287–0.438 — step prices make the last
price the right anchor); jump-boosted sigma after a big recent move
(worse, 0.277–0.287); Student-t quantile shapes (worse, 0.273–0.275);
a hard absolute leader threshold (price > 0.05) for the decay shift
(a touch better on the practice, 0.2722, but it would switch the shift
off entirely for the FAA family, whose prices sit below 0.05 — so the
scale-free ramp is used instead); per-group sigma multipliers
(0.2721, within noise of the chosen 0.2746, dropped as extra knobs
for no robust gain).

## Practice results (local scoring)

Scored against the next case's prices (p05 has no later case, so it is
not in the local numbers), 5 cases, 124,800 items, pinball loss in
percentage points:

| version | pooled | per-case mean | p00 | p01 | p02 | p03 | p04 |
|---------|--------|---------------|-----|-----|-----|-----|-----|
| house baseline | — | 0.398458 | 0.5361 | 0.3804 | 0.4084 | 0.3686 | 0.2988 |
| v2 (previous entry) | 0.318224 | 0.327831 | 0.5059 | 0.3092 | 0.3865 | 0.2366 | 0.2010 |
| v3 | 0.286033 | 0.291501 | 0.3534 | 0.3084 | 0.3314 | 0.2670 | 0.1974 |
| **v4 (this entry)** | **0.274599** | **0.279187** | **0.3120** | **0.3013** | **0.3201** | **0.2770** | **0.1855** |

v4 improves on v3 in 4 of 5 cases (p03 is slightly worse, 0.2770 vs
0.2670 — p03's horizon contains a +0.144 rise in a small market, which
no downward shift can anticipate). The pooled target of < 0.27 was NOT
reached (0.2746): after ~60 measured variants the score plateaus at
0.272–0.276, and the remaining loss is dominated by the median —
with an oracle median (same spreads) the pooled loss would be ~0.09,
so the gap is unpredictable news-driven moves of the leader markets,
not spread miscalibration.

**Holdout checks for v4's shift** (the shift coefficient was the one
parameter tuned on subsets): tuned on p00–p03, held-out p04 = 0.1855
(v3: 0.1974); tuned on p01–p04, held-out p00 = 0.3120 (v3: 0.3534).
Both held-out cases confirm the gain, so it is not a one-case artefact.

The house baseline's published practice figure is 0.363589 pooled over
all six cases. The previous entry's submitted figure, 0.279746, could
not be reproduced from the local data (its p00–p04 pooled is 0.318224,
so the difference must come from p05, which has no local reference
answers); v4 beats v3 by 0.011 pooled on the five verifiable cases, and
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
  ≤ 1 MiB (198–266 KB).
- A synthetic FAA-like case (4 markets, prices summing to 0.10, born at
  ~0.49 with a fast crash like the real FAA markets) was run to check
  the adaptive sum target and the scale-free shift: medians stay near
  the checkpoint prices (sum of medians at h240 = 0.0935 vs checkpoint
  sum 0.10), i.e. no pull toward a sum of 1, and the shift fires for
  the larger FAA-like markets (0.045) while leaving the smallest (0.010)
  untouched.
- Python 3.12.3 locally (grader uses 3.13); stdlib only, no network,
  no reads outside the case file.
