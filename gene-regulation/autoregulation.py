"""Self-activating gene: phase plane, fixed points and bistability.

Model (nondimensional rates, forward Euler):
    dr/dt = mu * p^h / (K^h + p^h) - chi_r * r      (mRNA)
    dp/dt = omega * r - chi_p * p                    (protein)

Reproduce everything in README.md with:  python autoregulation.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

MU, OMEGA, CHI_R, CHI_P, H = 1.0, 1.0, 1.0, 1.0, 2
K = 0.33
DT, T_END = 0.01, 20.0

FIG = Path(__file__).parent / "figures"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def hill(p, k=K, h=H):
    return p**h / (k**h + p**h)


def simulate(p0, r0, k=K, dt=DT, t_end=T_END):
    """Forward Euler; p0, r0 may be arrays (one trajectory per entry)."""
    p, r = np.array(p0, dtype=float), np.array(r0, dtype=float)
    n = int(round(t_end / dt))
    P, R = np.empty((n + 1,) + p.shape), np.empty((n + 1,) + p.shape)
    P[0], R[0] = p, r
    for i in range(n):
        dr = MU * hill(p, k) - CHI_R * r
        dp = OMEGA * r - CHI_P * p
        p, r = p + dt * dp, r + dt * dr
        P[i + 1], R[i + 1] = p, r
    return np.arange(n + 1) * dt, P, R


def fixed_points(k=K):
    """Steady states with h = 2 and unit rates: p = p^2 / (k^2 + p^2), r = p."""
    disc = 1 - 4 * k**2
    if disc < 0:
        return [0.0]
    return [0.0, (1 - np.sqrt(disc)) / 2, (1 + np.sqrt(disc)) / 2]


def is_stable(p, k=K):
    dhill = MU * H * k**H * p**(H - 1) / (k**H + p**H) ** 2
    J = np.array([[-CHI_R, dhill], [OMEGA, -CHI_P]])  # d(dr, dp)/d(r, p)
    return np.all(np.linalg.eigvals(J).real < 0)


def style(ax):
    ax.grid(color=GRID, lw=0.8)
    ax.tick_params(colors=MUTED)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    fps = fixed_points()
    for p in fps:
        print(f"fixed point p = r = {p:.4f}: {'stable' if is_stable(p) else 'unstable (saddle)'}")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    # (a) Hill function: cooperativity sharpens the switch
    s = np.linspace(0, 1.5, 300)
    for h, color in zip([1, 2, 10], [BLUE, ORANGE, AQUA]):
        axes[0].plot(s, hill(s, K, h), color=color, lw=2, label=f"h = {h}")
    axes[0].axvline(K, color=MUTED, ls=":", lw=1)
    axes[0].set_xlabel("activator (protein) concentration", color=MUTED)
    axes[0].set_ylabel("relative transcription rate", color=MUTED)
    axes[0].set_title("(a) Hill activation, K = 0.33", color=INK, fontsize=11)
    axes[0].legend(frameon=False)

    # (b) phase plane: trajectories coloured by the state they reach
    coarse = [a.ravel() for a in np.meshgrid(np.arange(0, 1.4, 0.1), np.arange(0, 1.4, 0.1))]
    fine = [a.ravel() for a in np.meshgrid(np.arange(0, 0.42, 0.03), np.arange(0, 0.42, 0.03))]
    p0, r0 = np.concatenate([coarse[0], fine[0]]), np.concatenate([coarse[1], fine[1]])
    t, P, R = simulate(p0, r0)
    high = P[-1] > fps[1]
    ax = axes[1]
    for k in range(len(p0)):
        ax.plot(P[:, k], R[:, k], color=ORANGE if high[k] else BLUE, lw=0.8, alpha=0.7)
    pp = np.linspace(0, 1.4, 300)
    ax.plot(pp, pp, color=INK, lw=1.5, ls="--", label="dp/dt = 0  (r = p)")
    ax.plot(pp, hill(pp), color=INK, lw=1.5, label="dr/dt = 0  (r = Hill(p))")
    for p in fps:
        stable = is_stable(p)
        ax.plot(p, p, "o", ms=9, mec=INK, mfc=INK if stable else "white", zorder=5)
    ax.plot([], [], color=ORANGE, lw=1.5, label="→ high-expression state")
    ax.plot([], [], color=BLUE, lw=1.5, label="→ off state")
    ax.set_xlim(0, 1.35)
    ax.set_ylim(0, 1.35)
    ax.set_xlabel("protein p", color=MUTED)
    ax.set_ylabel("mRNA r", color=MUTED)
    ax.set_title("(b) phase plane: two stable states", color=INK, fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    # (c) the switch: number of steady states vs K
    ks = np.linspace(0.05, 0.7, 400)
    for k in ks:
        for p in fixed_points(k):
            axes[2].plot(k, p, ".", ms=2.5, color=INK if is_stable(p, k) else ORANGE)
    axes[2].axvline(K, color=MUTED, ls=":", lw=1)
    axes[2].axvline(0.5, color=MUTED, ls="--", lw=1)
    axes[2].text(0.51, 0.85, "saddle-node\nK = 0.5", color=MUTED, fontsize=9)
    axes[2].plot([], [], ".", color=INK, label="stable")
    axes[2].plot([], [], ".", color=ORANGE, label="unstable")
    axes[2].set_xlabel("activation threshold K", color=MUTED)
    axes[2].set_ylabel("steady-state protein p*", color=MUTED)
    axes[2].set_title("(c) bistability is lost for K > 0.5", color=INK, fontsize=11)
    axes[2].legend(frameon=False, fontsize=9)

    for ax in axes:
        style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "autoregulation.png", dpi=200)

    print("\nfinal state from the four initial conditions of the original exercise:")
    for p_init, r_init in [(0.4, 0.0), (0.0, 0.4), (0.2, 0.0), (0.0, 0.2), (0.01, 0.01), (0.4, 0.4)]:
        _, P1, R1 = simulate(p_init, r_init)
        print(f"  (p, r) = ({p_init}, {r_init}) -> p = {P1[-1]:.3f}")
