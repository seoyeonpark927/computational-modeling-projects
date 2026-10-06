"""Perceptron learning: Boolean functions and the storage capacity of random patterns.

Rule (one pattern at a time, random order):  w <- w + eta * (y_target - sign(w.x)) * x
Reproduce everything in README.md with:  python perceptron.py      (~1 min)
"""

from math import comb
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FIG = Path(__file__).parent / "figures"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

# inputs (x1, x2, bias = -1) and targets for the Boolean functions
X_BOOL = np.array([[1, 1, -1], [1, -1, -1], [-1, 1, -1], [-1, -1, -1]])
TARGETS = {"AND": [1, -1, -1, -1], "OR": [1, 1, 1, -1], "XOR": [-1, 1, 1, -1]}


def train_online(X, y, n_presentations, rng, eta=1.0):
    """Return the 0/1 correctness of each presentation (before the update)."""
    w = np.zeros(X.shape[1])
    correct = np.empty(n_presentations, dtype=int)
    for i, k in enumerate(rng.integers(0, len(X), n_presentations)):
        out = np.sign(X[k] @ w)
        correct[i] = out == y[k]
        if not correct[i]:
            w += eta * (y[k] - out) * X[k]
    return correct


def learns_all(X, y, max_epochs, rng, eta=1.0):
    """Train in random-order epochs; True once every pattern is classified correctly."""
    w = np.zeros(X.shape[1])
    for _ in range(max_epochs):
        for k in rng.permutation(len(X)):
            out = np.sign(X[k] @ w)
            if out != y[k]:
                w += eta * (y[k] - out) * X[k]
        if np.all(np.sign(X @ w) == y):
            return True
    return False


def cover_probability(M, N):
    """Fraction of random dichotomies of M points in general position in N dims
    that are linearly separable through the origin (Cover, 1965)."""
    return 2.0 ** (1 - M) * sum(comb(M - 1, k) for k in range(min(N, M)))


def random_problem(M, N, rng):
    X = rng.choice([-1, 1], size=(M, N))
    y = rng.permutation(np.repeat([1, -1], [M - M // 2, M // 2]))
    return X, y


def style(ax):
    ax.grid(color=GRID, lw=0.8)
    ax.tick_params(colors=MUTED)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    rng = np.random.default_rng(1)

    # (a) Boolean functions: moving accuracy over 500 presentations, 50 runs
    window = 20
    curves = {}
    for name, y in TARGETS.items():
        runs = np.array([train_online(X_BOOL, np.array(y), 500, rng) for _ in range(50)])
        curves[name] = runs.mean(0)
        print(f"{name}: accuracy over the last 100 presentations = {runs[:, -100:].mean():.2f}")

    # (b) capacity: fraction of random problems learned perfectly, N = 50
    N, runs, max_epochs = 50, 20, 500
    Ms = np.arange(10, 161, 10)
    learned = []
    for M in Ms:
        ok = np.mean([learns_all(*random_problem(M, N, rng), max_epochs, rng) for _ in range(runs)])
        learned.append(ok)
        print(f"M = {M:3d} (alpha = {M / N:.1f}): learned {ok:.2f}, Cover {cover_probability(M, N):.2f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    ax = axes[0]
    for (name, curve), color, ls in zip(curves.items(), [BLUE, AQUA, ORANGE], ["--", "-", "-"]):
        smooth = np.convolve(curve, np.ones(window) / window, mode="valid")
        ax.plot(np.arange(len(smooth)) + window, smooth, color=color, lw=2, ls=ls, zorder=3 if ls == "--" else 2,
                label=name + (" (overlaps OR)" if name == "AND" else ""))
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("pattern presentations", color=MUTED)
    ax.set_ylabel("fraction correct (50 runs)", color=MUTED)
    ax.set_title("(a) XOR is not linearly separable", color=INK, fontsize=11)
    ax.legend(frameon=False)

    ax = axes[1]
    alpha = np.linspace(0.2, 3.2, 300)
    ax.plot(alpha, [cover_probability(int(round(a * N)), N) for a in alpha],
            color=MUTED, lw=2, ls="--", label="Cover's theorem")
    ax.plot(Ms / N, learned, "o", color=BLUE, ms=7, label=f"perceptron, {max_epochs} epochs, {runs} runs per point")
    ax.axvline(2, color=MUTED, lw=1, ls=":")
    ax.set_xlabel("load α = patterns / inputs  (N = 50)", color=MUTED)
    ax.set_ylabel("fraction of problems learned", color=MUTED)
    ax.set_title("(b) capacity: random patterns are separable up to α ≈ 2", color=INK, fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="center left")
    for ax in axes:
        style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "perceptron.png", dpi=200)
