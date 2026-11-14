#!/usr/bin/env python3
"""Forecast Polymarket Yes prices 240 minutes ahead, five quantiles per minute.

Method (v3): a random-walk forecast around each market's checkpoint price,
with the spread fitted to the market's own recent volatility and the median
pulled toward the event's own recent probability-mass level.

* Median: the checkpoint price, plus a drift that pulls the whole event's
  prices back toward the level their sum sat at over the last 12 hours.
  The markets of one event are alternative answers to the same question,
  so their Yes prices together are one probability mass: when the sum at
  the checkpoint sits above (or below) its recent level, the excess is
  traded away over the next hour, distributed across markets in proportion
  to each market's share of the sum. The target is the MEDIAN of the sum
  over the last 720 minutes -- a level the event itself set, not a number
  assumed in advance -- so the same rule works whether the family's prices
  add up to about 1 (a proper mutually exclusive event) or to about 0.1
  (a family of long-shot attributions). The median, not the mean, is used
  because a market's first minutes trade near 0.5 before finding its
  level, which spikes the sum briefly; the median ignores those spikes.
  Strength 0.8, time constant 60 minutes.
* Spread: each market's own per-minute volatility, the LARGER of a fast
  (half-life 30 minutes) and a slow (half-life 2 hours) exponentially
  weighted variance of its one-minute price changes (last 24 hours,
  winsorised at +-0.05), scaled by 0.3 and floored at 0.0002. Volatility
  clusters: it jumps up fast after a burst and decays slowly, so the
  spread should rise fast and fall slow. The spread grows like sqrt(h)
  for minute h, with the 5%/95% quantiles widened a further 1.3x because
  the moves are heavy-tailed (quantiles of the standard normal).
* No drift from a market's own recent trend: one-minute moves show no
  usable autocorrelation, and pulling the median toward a trailing mean
  made the practice score worse (old deviations are regime shifts, not
  temporary excursions).

Quantiles are clamped to [0, 1], rounded to 6 decimals and forced
non-decreasing. Standard library only, no network, deterministic.

Practice loss (local scoring, pinball loss in percentage points):
p00-p04, scored against the next case's prices -- old entry 0.318224
pooled, new entry 0.286033 pooled (per-case 0.3534 / 0.3084 / 0.3314 /
0.2670 / 0.1974, mean 0.2915), versus the old entry's 0.5059 / 0.3092 /
0.3865 / 0.2366 / 0.2010 (mean 0.3278). The house baseline scores
0.363589 pooled over all six practice cases. Case p05 has no later case
to score it against locally, so it is not in the local numbers above;
the new entry's p05 drift is near zero (its sum sits almost exactly on
its recent level), and on a frozen-price approximation of p05 it ties
the old entry, so the six-case total is expected to improve by about
the same margin. The practice is a different event, so this is a guide
only, not a promise about the final.
"""
import json
import math
import os
from datetime import datetime, timezone

INPUT_DIR = os.environ.get("MAKEFUN_INPUT_DIR", "/input")
OUTPUT_DIR = os.environ.get("MAKEFUN_OUTPUT_DIR", "/output")

QUANTILES = [0.05, 0.25, 0.5, 0.75, 0.95]
# The standard normal's quantile at each level.
Z = [-1.644854, -0.67449, 0.0, 0.67449, 1.644854]

# Half-lives, in minutes, of the fast and slow exponential weights on
# past one-minute moves. The max of the two variances is the volatility.
SIGMA_HALF_LIFE_FAST_MIN = 30.0
SIGMA_HALF_LIFE_SLOW_MIN = 120.0
# One-minute moves are winsorised to this range before entering the variance.
DIFF_CAP = 0.05
# Only the most recent that many minutes of history feed the volatility.
SIGMA_WINDOW_MIN = 1440
# The raw per-minute volatility overstates the horizon's moves, so scale it.
SIGMA_MULT = 0.3
# No market is ever forecast as perfectly still.
SIGMA_FLOOR = 0.0002
# Extra widening of the outer (5% / 95%) quantiles: heavy-tailed moves.
TAIL_MULT = 1.3
# Fraction of the sum's excess over its recent level removed over the horizon.
SUM_REVERT = 0.8
# Time constant, in minutes, of the sum reversion.
SUM_TAU_MIN = 60.0
# The sum's recent level is the median of its last that many minutes.
SUM_TARGET_WINDOW_MIN = 720


def parse_ts(s):
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%S.%fZ").replace(
        tzinfo=timezone.utc)


def price6(x):
    """x clamped to [0, 1] and rounded to 6 decimals (0, never -0)."""
    return round(min(1.0, max(0.0, x)), 6) + 0.0


def median(values):
    s = sorted(values)
    n = len(s)
    if n == 0:
        return 0.0
    if n % 2:
        return s[n // 2]
    return 0.5 * (s[n // 2 - 1] + s[n // 2])


def sum_target(case):
    """Median of the markets' summed Yes price over the last 12 hours.

    The sum is built on the absolute minute grid: a market contributes
    from its first minute. The median over the window is robust to the
    brief spike a newly created market causes (it first trades near 0.5
    before finding its level) and to markets joining mid-window.
    Returns None when there is too little history to judge a level.
    """
    markets = case["markets"]
    if not markets:
        return None
    t0 = min(parse_ts(m["start"]) for m in markets)
    t1 = parse_ts(case["at"])
    nmin = int((t1 - t0).total_seconds() // 60) + 1
    if nmin < 60:
        return None
    sums = [0.0] * nmin
    for m in markets:
        off = int((parse_ts(m["start"]) - t0).total_seconds() // 60)
        prices = m["prices"]
        for i, v in enumerate(prices):
            sums[off + i] += v
    window = sums[-SUM_TARGET_WINDOW_MIN:]
    return median(window)


def market_sigma(prices):
    """Per-minute volatility: the larger of a fast and a slow EWMA variance.

    Rises fast after a burst (the 30-minute half-life reacts within
    minutes) and falls slow (the 2-hour half-life holds the regime level).
    Both variances are computed in a single pass over the window, then
    scaled and floored.
    """
    n = len(prices)
    if n < 2:
        return SIGMA_FLOOR
    alpha_fast = 1.0 - 0.5 ** (1.0 / SIGMA_HALF_LIFE_FAST_MIN)
    alpha_slow = 1.0 - 0.5 ** (1.0 / SIGMA_HALF_LIFE_SLOW_MIN)
    var_fast = 0.0
    var_slow = 0.0
    start = n - SIGMA_WINDOW_MIN if n > SIGMA_WINDOW_MIN else 1
    for i in range(start, n):
        d = prices[i] - prices[i - 1]
        if d > DIFF_CAP:
            d = DIFF_CAP
        elif d < -DIFF_CAP:
            d = -DIFF_CAP
        d2 = d * d
        var_fast = alpha_fast * d2 + (1.0 - alpha_fast) * var_fast
        var_slow = alpha_slow * d2 + (1.0 - alpha_slow) * var_slow
    var = var_fast if var_fast > var_slow else var_slow
    sigma = math.sqrt(var) * SIGMA_MULT
    return sigma if sigma > SIGMA_FLOOR else SIGMA_FLOOR


def forecast_case(case):
    horizon = int(case["horizon_minutes"])
    markets = case["markets"]
    total = 0.0
    for m in markets:
        total += m["prices"][-1] if m["prices"] else 0.0
    target = sum_target(case)
    excess = (total - target) if target is not None else 0.0
    sqrt_h = [math.sqrt(h + 1) for h in range(horizon)]
    revert_h = [1.0 - math.exp(-(h + 1) / SUM_TAU_MIN) for h in range(horizon)]
    forecasts = []
    for m in markets:
        prices = m["prices"]
        if not prices:
            # No history at all: stay at 0.5 with the floor spread.
            prices = [0.5]
        x = prices[-1]
        sigma = market_sigma(prices)
        share = x / total if total > 0.0 else 0.0
        drift_total = -SUM_REVERT * excess * share
        rows = []
        for h in range(horizon):
            med = x + drift_total * revert_h[h]
            spread = sigma * sqrt_h[h]
            row = []
            for j, z in enumerate(Z):
                zz = z * TAIL_MULT if j == 0 or j == 4 else z
                row.append(price6(med + zz * spread))
            for j in range(1, 5):
                if row[j] < row[j - 1]:
                    row[j] = row[j - 1]
            rows.append(row)
        forecasts.append({"id": m["id"], "q": rows})
    return forecasts


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    inputs = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith(".in"))
    for name in inputs:
        case_id = name[: -len(".in")]
        with open(os.path.join(INPUT_DIR, name), "r", encoding="utf-8") as f:
            case = json.load(f)
        if len(case["quantiles"]) != len(QUANTILES) or any(
            t != q for t, q in zip(case["quantiles"], QUANTILES)
        ):
            raise ValueError(f"{name} asks for quantiles {case['quantiles']}, not {QUANTILES}")
        forecasts = forecast_case(case)
        out = json.dumps({"forecasts": forecasts}, separators=(",", ":")) + "\n"
        with open(os.path.join(OUTPUT_DIR, case_id + ".out"), "w", encoding="utf-8") as f:
            f.write(out)


if __name__ == "__main__":
    main()
