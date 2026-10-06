"""Hodgkin-Huxley and leaky integrate-and-fire neurons, forward Euler.

Units: time in ms, voltage in mV.
  HH : conductances in nS, capacitance in pF, current in pA  (nS*mV = pA, pA/pF = mV/ms)
  LIF: conductance in uS, capacitance in nF, current in nA  (uS*mV = nA, nA/nF = mV/ms)

All functions accept an array of input currents and simulate one neuron per
current in parallel.
"""

import numpy as np

# ---------------------------------------------------------------- Hodgkin-Huxley
HH = dict(G_Na=400.0, G_K=200.0, G_L=2.0, E_Na=99.0, E_K=-85.0, V_L=-65.0, C=2.0)
HH_SPIKE_THRESHOLD = -20.0  # mV, upward crossing counts as a spike


def rates(V):
    """Opening (alpha) and closing (beta) rates of the m, h, n gates, in 1/ms."""
    a_m = 0.1 * (V + 40) / (1 - np.exp(-0.1 * (V + 40)))
    b_m = 4 * np.exp(-0.0556 * (V + 65))
    a_h = 0.07 * np.exp(-0.05 * (V + 65))
    b_h = 1 / (1 + np.exp(-0.1 * (V + 35)))
    a_n = 0.01 * (V + 55) / (1 - np.exp(-0.1 * (V + 55)))
    b_n = 0.125 * np.exp(-0.0125 * (V + 65))
    return a_m, b_m, a_h, b_h, a_n, b_n


def simulate_hh(I0, dt=0.01, t_max=200.0, t_on=40.0, record=False):
    """Step current I0 (pA, scalar or array) switched on at t_on.

    Returns spike times per neuron, and full traces if record=True.
    """
    p = HH
    I0 = np.atleast_1d(np.asarray(I0, dtype=float))
    t = np.arange(0, t_max + dt / 2, dt)
    V = np.full(I0.shape, p["V_L"])
    a_m, b_m, a_h, b_h, a_n, b_n = rates(V)
    m, h, n = a_m / (a_m + b_m), a_h / (a_h + b_h), a_n / (a_n + b_n)

    spikes = [[] for _ in I0]
    trace = {k: np.empty((len(t), len(I0))) for k in ["V", "g_Na", "g_K", "I_m", "I_e"]} if record else None

    with np.errstate(all="ignore"):  # large dt makes the integration blow up on purpose
        for i, ti in enumerate(t):
            I_e = np.where(ti < t_on, 0.0, I0)
            g_Na = p["G_Na"] * m**3 * h
            g_K = p["G_K"] * n**4
            I_m = g_Na * (V - p["E_Na"]) + g_K * (V - p["E_K"]) + p["G_L"] * (V - p["V_L"])
            if record:
                for k, val in zip(trace, [V, g_Na, g_K, I_m, I_e]):
                    trace[k][i] = val

            a_m, b_m, a_h, b_h, a_n, b_n = rates(V)
            V_new = V + dt * (I_e - I_m) / p["C"]
            m = m + dt * (a_m * (1 - m) - b_m * m)
            h = h + dt * (a_h * (1 - h) - b_h * h)
            n = n + dt * (a_n * (1 - n) - b_n * n)

            crossed = (V < HH_SPIKE_THRESHOLD) & (V_new >= HH_SPIKE_THRESHOLD)
            for k in np.flatnonzero(crossed):
                spikes[k].append(ti + dt)
            V = V_new
    return (t, spikes, trace) if record else spikes


# --------------------------------------------------------- leaky integrate-and-fire
LIF = dict(C=1.0, G_L=0.05, V_L=-65.0, V_spk=-45.0, V_reset=-65.0, t_ref=2.0)


def simulate_lif(I0, dt=0.1, t_max=200.0, t_on=40.0, record=False):
    """Step current I0 (nA, scalar or array) switched on at t_on."""
    p = LIF
    I0 = np.atleast_1d(np.asarray(I0, dtype=float))
    t = np.arange(0, t_max + dt / 2, dt)
    V = np.full(I0.shape, p["V_L"])
    refractory = np.zeros(I0.shape, dtype=int)
    spikes = [[] for _ in I0]
    V_trace = np.empty((len(t), len(I0))) if record else None

    for i, ti in enumerate(t):
        if record:
            V_trace[i] = V
        I_e = np.where(ti < t_on, 0.0, I0)
        free = refractory == 0
        V = np.where(free, V + dt * (-p["G_L"] * (V - p["V_L"]) + I_e) / p["C"], p["V_reset"])
        refractory = np.maximum(refractory - 1, 0)
        fired = free & (V >= p["V_spk"])
        for k in np.flatnonzero(fired):
            spikes[k].append(ti + dt)
        V = np.where(fired, p["V_reset"], V)
        refractory = np.where(fired, int(round(p["t_ref"] / dt)), refractory)
        if record:
            V_trace[i] = np.where(fired, 20.0, V_trace[i])  # draw the spike
    return (t, spikes, V_trace) if record else spikes


# ------------------------------------------------------------------- analysis
def firing_rate(spike_times, t_start):
    """Sustained rate (Hz) from the mean inter-spike interval after t_start.

    Neurons that fire a few spikes and then stop get rate 0.
    """
    s = np.asarray(spike_times)
    s = s[s >= t_start]
    if len(s) < 3:
        return 0.0
    return 1000.0 / np.mean(np.diff(s))
