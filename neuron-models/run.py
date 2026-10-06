"""Reproduce every figure and number in neuron-models/README.md.

    python run.py          (~2 min)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from neurons import HH, LIF, firing_rate, simulate_hh, simulate_lif

FIG = Path(__file__).parent / "figures"
FIG.mkdir(exist_ok=True)
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def style(ax):
    ax.grid(color=GRID, lw=0.8)
    ax.tick_params(colors=MUTED)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)


# 1. One HH neuron driven by a 200 pA step ----------------------------------
t, spikes, tr = simulate_hh(200.0, dt=0.01, t_max=100.0, t_on=40.0, record=True)
fig, axes = plt.subplots(4, 1, figsize=(7, 6.5), sharex=True)
for ax, (key, label, color) in zip(axes, [
        ("V", "V (mV)", BLUE), ("g_Na", "g_Na (nS)", ORANGE),
        ("g_K", "g_K (nS)", AQUA), ("I_e", "I_e (pA)", INK)]):
    ax.plot(t, tr[key][:, 0], color=color, lw=1.5)
    ax.set_ylabel(label, color=MUTED)
    style(ax)
axes[-1].set_xlabel("time (ms)", color=MUTED)
axes[0].set_title("Hodgkin–Huxley neuron, 200 pA step at t = 40 ms", color=INK)
fig.tight_layout()
fig.savefig(FIG / "hh_spike_train.png", dpi=200)

# 2. Thresholds of the HH neuron (as in the original analysis: 200 ms, dt 0.01)
I_scan = np.arange(0.0, 200.0, 0.5)
counts = np.array([len(s) for s in simulate_hh(I_scan, dt=0.01, t_max=200.0)])
I_single = I_scan[np.argmax(counts >= 1)]
I_repetitive = I_scan[np.argmax(counts >= 6)]
print(f"HH: first action potential at I = {I_single} pA, "
      f"repetitive firing from I = {I_repetitive} pA")

# 3. f-I curves --------------------------------------------------------------
T_MAX, T_ON = 1000.0, 40.0
I_hh = np.arange(0.0, 405.0, 5.0)
rate_hh = [firing_rate(s, T_ON + 200) for s in simulate_hh(I_hh, dt=0.01, t_max=T_MAX)]

I_lif = np.linspace(0.0, 3.0, 61)
rate_lif = [firing_rate(s, T_ON + 200) for s in simulate_lif(I_lif, dt=0.1, t_max=T_MAX)]
tau = LIF["C"] / LIF["G_L"]
I_rh = LIF["G_L"] * (LIF["V_spk"] - LIF["V_L"])
I_fine = np.linspace(I_rh * 1.0001, 3.0, 400)
rate_theory = 1000.0 / (LIF["t_ref"] + tau * np.log(
    (I_fine / LIF["G_L"]) / (I_fine / LIF["G_L"] - (LIF["V_spk"] - LIF["V_L"]))))
print(f"LIF: theoretical rheobase {I_rh:.2f} nA, membrane time constant {tau:.0f} ms")

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
axes[0].plot(I_hh, rate_hh, "o", color=BLUE, ms=5)
axes[0].set_xlabel("input current (pA)", color=MUTED)
axes[0].set_title("Hodgkin–Huxley: discontinuous (type II)", color=INK, fontsize=11)
axes[1].plot(I_fine, rate_theory, color=MUTED, lw=2, ls="--", label="analytical")
axes[1].plot(I_lif, rate_lif, "o", color=ORANGE, ms=5, label="simulation (dt = 0.1 ms)")
axes[1].axvline(I_rh, color=MUTED, lw=1, ls=":")
axes[1].set_xlabel("input current (nA)", color=MUTED)
axes[1].set_title("Leaky integrate-and-fire: continuous (type I)", color=INK, fontsize=11)
axes[1].legend(frameon=False, fontsize=9)
for ax in axes:
    ax.set_ylabel("firing rate (Hz)", color=MUTED)
    style(ax)
fig.tight_layout()
fig.savefig(FIG / "fi_curves.png", dpi=200)

# 4. Time step and stability of forward Euler for HH -------------------------
fig, ax = plt.subplots(figsize=(7, 3.6))
for dt, color in zip([0.01, 0.025, 0.03], [BLUE, ORANGE, AQUA]):
    t, spk, tr = simulate_hh(200.0, dt=dt, t_max=100.0, record=True)
    V = tr["V"][:, 0]
    finite = np.isfinite(V) & (np.abs(V) < 500)
    end = np.argmin(finite) if not finite.all() else len(V)
    note = f" ({len(spk[0])} spikes)" if finite.all() else f" (diverges at t = {t[end]:.1f} ms)"
    ax.plot(t[:end], V[:end], color=color, lw=1.5, label=f"dt = {dt} ms{note}")
    print(f"HH dt={dt}: {'stable' if finite.all() else 'diverges'}")
ax.set_ylim(-100, 80)
ax.set_xlabel("time (ms)", color=MUTED)
ax.set_ylabel("V (mV)", color=MUTED)
ax.set_title("Forward Euler time step and the Hodgkin–Huxley model", color=INK)
ax.legend(frameon=False, fontsize=9, loc="lower right")
style(ax)
fig.tight_layout()
fig.savefig(FIG / "euler_stability.png", dpi=200)
print("figures saved to", FIG)
