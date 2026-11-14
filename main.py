#!/usr/bin/env python3
"""Forecast Polymarket Yes prices 240 minutes ahead, five quantiles per minute.

Method (v2): a random-walk forecast around the checkpoint price, with the
spread fitted to each market's own recent behaviour.

* Spread: the market's recent per-minute volatility, as the exponentially
  weighted root-mean-square of its one-minute price changes (half-life
  30 minutes, last 12 hours), floored at 0.0008 and capped at 0.03 price
  units per square-root minute. Recent moves count more than old ones, so
  the spread tracks the market's current liveliness.
* Shape: quantiles of a Student-t distribution with 4 degrees of freedom
  (z = -2.131847, -0.740697, 0, 0.740697, 2.131847), scaled by 0.4 and by
  sqrt(h) for minute h. Prices move in steps -- mostly flat, then a tick,
  now and then a big jump -- so the forecast is narrow in the middle (cheap
  on quiet minutes) with heavy tails (cheap when a jump lands outside it).
  The 0.4 scale was chosen so the 90% interval covers about 90% of minutes
  on the practice cases.
* Median: the checkpoint price itself. Short-term drift was tested on the
  practice cases and made the median strictly worse in every form tried
  (mean-reverting step prices punish chasing the last move), so no drift
  is applied.

Quantiles are clamped to [0, 1], rounded to 6 decimals and forced
non-decreasing. Standard library only, no network, deterministic.

Practice loss (6 cases, scored locally against the next case's prices,
p05 against Polymarket's public history): 0.2797, versus the house
baseline's 0.3636. The practice is a different event, so this is a guide
only, not a promise about the final.
"""
import json
import math
import os

INPUT_DIR = os.environ.get("MAKEFUN_INPUT_DIR", "/input")
OUTPUT_DIR = os.environ.get("MAKEFUN_OUTPUT_DIR", "/output")

QUANTILES = [0.05, 0.25, 0.5, 0.75, 0.95]
# Student-t(4 df) quantiles at the five levels.
Z = [-2.131847, -0.740697, 0.0, 0.740697, 2.131847]

SIGMA_FLOOR = 0.0008
SIGMA_CAP = 0.03
SIGMA_MULT = 0.4
VOL_HALF_LIFE = 30.0
VOL_WINDOW = 720


def price6(x):
    v = round(min(1.0, max(0.0, x)) * 1e6) / 1e6
    return v + 0.0  # avoid -0.0


def recent_sigma(prices):
    """EWMA root-mean-square of one-minute price changes."""
    n = len(prices)
    if n < 2:
        return SIGMA_FLOOR
    diffs = [prices[i + 1] - prices[i] for i in range(n - 1)]
    a = 1.0 - 0.5 ** (1.0 / VOL_HALF_LIFE)
    var = 0.0
    wsum = 0.0
    w = 1.0
    for d in reversed(diffs[-VOL_WINDOW:]):
        var += w * d * d
        wsum += w
        w *= (1.0 - a)
    var /= wsum if wsum else 1.0
    sigma = math.sqrt(var)
    return min(SIGMA_CAP, max(SIGMA_FLOOR, sigma))


def forecast_market(market, horizon):
    prices = market["prices"]
    if not prices:
        raise ValueError("market %s has no prices" % market["id"])
    x = prices[-1]
    sigma = recent_sigma(prices) * SIGMA_MULT
    rows = []
    for h in range(1, horizon + 1):
        s = sigma * math.sqrt(h)
        row = [price6(x + z * s) for z in Z]
        # enforce non-decreasing after rounding
        for j in range(1, 5):
            if row[j] < row[j - 1]:
                row[j] = row[j - 1]
        rows.append(row)
    return rows


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    files = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith(".in"))
    for fname in files:
        case_id = fname[:-3]
        with open(os.path.join(INPUT_DIR, fname), "r", encoding="utf-8") as f:
            case = json.load(f)
        horizon = int(case["horizon_minutes"])
        forecasts = [
            {"id": m["id"], "q": forecast_market(m, horizon)}
            for m in case["markets"]
        ]
        out_path = os.path.join(OUTPUT_DIR, case_id + ".out")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(json.dumps({"forecasts": forecasts}) + "\n")


if __name__ == "__main__":
    main()
