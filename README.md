# FAA outage attribution — 4-hour Polymarket price forecast

Entry for the make.fun competition `polymarket-faa-attribution-02`.
`python3 main.py` reads every `*.in` case in `$MAKEFUN_INPUT_DIR` (or
`/input`) and writes one `*.out` forecast per case to `$MAKEFUN_OUTPUT_DIR`
(or `/output`). Standard library only, no network, deterministic.

## Approach

A random-walk forecast around each market's checkpoint price:

- **Spread fitted per market.** The per-minute volatility is the
  exponentially weighted root-mean-square of the market's one-minute price
  changes (half-life 30 minutes, last 12 hours), floored at 0.0008 and
  capped at 0.03 price units per square-root minute. Recent moves count
  more, so the spread tracks how lively the market is right now.
- **Heavy-tailed shape.** Quantiles follow a Student-t distribution with
  4 degrees of freedom (z = -2.131847, -0.740697, 0, 0.740697, 2.131847),
  scaled by 0.4 and by sqrt(h) for minute h. Prices move in steps — mostly
  flat, then a tick, now and then a big jump — so the forecast is narrow in
  the middle (cheap on quiet minutes) with heavy tails (cheap when a jump
  lands outside it). The 0.4 scale makes the 90% interval cover about 90%
  of minutes on the practice cases.
- **Median = checkpoint price, no drift.** Short-term drift was tested on
  the practice cases and made the median strictly worse in every form
  tried (step prices with mean reversion punish chasing the last move), so
  no drift is applied.

Quantiles are clamped to [0, 1], rounded to 6 decimals, and forced
non-decreasing.

## Practice results (local scoring)

Scored against the next case's prices (p05 against Polymarket's public
price history), 6 cases, 153,600 items:

| version | loss | median error | 90% coverage |
|---------|------|--------------|--------------|
| house baseline | 0.363589 | 0.767503 | 0.949642 |
| v1 (EWMA sigma + drift) | 0.424386 | 1.098500 | 0.907487 |
| **v2 (this entry)** | **0.279746** | **0.767503** | **0.911230** |

The practice is cut from a different event (Next CEO of NVIDIA), so these
numbers are a guide only, not a promise about the final.
