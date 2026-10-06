"""Poisson spike trains: raster, PSTH, ISI variability (CV) and count variability (Fano factor).

Spike trains are built from exponential inter-spike intervals, ISI = -ln(u) / rate.
Reproduce everything in README.md with:  python spike_trains.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RATE, N_TRIALS, T = 10.0, 50, 5.0       # spikes/s, trials, seconds
FIG = Path(__file__).parent / "figures"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
BLUE, ORANGE = "#2a78d6", "#eb6834"


def poisson_trains(rate, n_trials, duration, rng):
    """List of spike-time arrays, one per trial."""
    n_isi = int(2 * rate * duration) + 20  # enough intervals to cover the duration
    trains = []
    for _ in range(n_trials):
        t = np.cumsum(-np.log(rng.random(n_isi)) / rate)
        while t[-1] < duration:  # extremely rare: extend if still short
            t = np.concatenate([t, t[-1] + np.cumsum(-np.log(rng.random(n_isi)) / rate)])
        trains.append(t[t <= duration])
    return trains


def psth(trains, duration, bin_width):
    edges = np.arange(0, duration + bin_width / 2, bin_width)
    counts = sum(np.histogram(tr, edges)[0] for tr in trains)
    return (edges[:-1] + edges[1:]) / 2, counts / (len(trains) * bin_width)


def cv(train):
    isi = np.diff(train)
    return np.std(isi) / np.mean(isi) if len(isi) > 1 else np.nan


def fano(trains):
    counts = np.array([len(tr) for tr in trains])
    return np.var(counts) / np.mean(counts)


def style(ax):
    ax.grid(color=GRID, lw=0.8)
    ax.tick_params(colors=MUTED)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    rng = np.random.default_rng(0)
    trains = poisson_trains(RATE, N_TRIALS, T, rng)

    rate_hat = sum(len(tr) for tr in trains) / (N_TRIALS * T)
    cvs = np.array([cv(tr) for tr in trains])
    print(f"mean rate {rate_hat:.2f} Hz (theory {RATE}), mean CV {np.nanmean(cvs):.3f} (theory 1), "
          f"Fano factor {fano(trains):.3f} (theory 1)")

    # Variability of the estimates vs. recording length, 20 repeats per T
    durations = [10, 100, 200, 400, 800, 1600, 3200, 6400]
    repeats = 20
    cv_T = np.array([[cv(poisson_trains(RATE, 1, d, rng)[0]) for _ in range(repeats)] for d in durations])
    ff_T = np.array([[fano(poisson_trains(RATE, N_TRIALS, d, rng)) for _ in range(repeats)] for d in durations])
    for d, c, f in zip(durations, cv_T, ff_T):
        print(f"T = {d:5d} s: CV {c.mean():.3f} ± {c.std():.3f}, FF {f.mean():.3f} ± {f.std():.3f}")

    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    ax = axes[0]
    for k, tr in enumerate(trains):
        ax.plot(tr, np.full(len(tr), k + 1), "|", color=INK, ms=3)
    ax.set_xlabel("time (s)", color=MUTED)
    ax.set_ylabel("trial", color=MUTED)
    ax.set_title(f"(a) {N_TRIALS} Poisson spike trains, {RATE:g} Hz", color=INK, fontsize=11)

    centers, f_bin = psth(trains, T, 0.2)
    axes[1].step(centers, f_bin, where="mid", color=BLUE, lw=1.8)
    axes[1].axhline(RATE, color=MUTED, ls="--", lw=1.2, label="true rate")
    axes[1].set_ylim(0, 2 * RATE)
    axes[1].set_xlabel("time (s)", color=MUTED)
    axes[1].set_ylabel("firing rate (Hz)", color=MUTED)
    axes[1].set_title("(b) PSTH, 200 ms bins", color=INK, fontsize=11)
    axes[1].legend(frameon=False, fontsize=9)

    ax = axes[2]
    for data, color, label in [(cv_T, BLUE, "CV of ISIs (one train)"),
                               (ff_T, ORANGE, f"Fano factor ({N_TRIALS} trials)")]:
        ax.errorbar(durations, data.mean(1), yerr=data.std(1), fmt="o-", color=color,
                    capsize=3, lw=1.8, ms=6, label=label)
    ax.axhline(1, color=MUTED, ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_ylim(0.5, 1.5)
    ax.set_xlabel("recording length T (s)", color=MUTED)
    ax.set_ylabel(f"value (mean ± SD, {repeats} repeats)", color=MUTED)
    ax.set_title("(c) CV converges with T; Fano factor does not", color=INK, fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    for ax in axes:
        style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "spike_train_statistics.png", dpi=200)
