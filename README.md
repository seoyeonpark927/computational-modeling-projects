# Computational modeling projects

Short, self-contained projects in computational neuroscience, systems biology and
numerical analysis. Each one asks a single question, answers it with a simulation,
and checks the result against theory where possible.

| Project | Question | Key result |
|---|---|---|
| [Spiking neuron models](neuron-models/) | How do Hodgkin–Huxley and integrate-and-fire neurons convert current into firing rate? | HH jumps to ~75 Hz at onset (type II), LIF rises continuously and matches its analytical f–I curve. Forward Euler for HH is only usable for dt ≲ 0.02 ms. |
| [Gene regulation](gene-regulation/) | Can a self-activating gene act as a memory? | Yes: three steady states, two stable. Bistability is lost in a saddle-node bifurcation at K = 0.5. |
| [Spike-train statistics](spike-train-statistics/) | How much data does it take to measure neural variability? | The ISI CV converges with recording length. The Fano factor's precision is set by the number of trials, not by duration. |
| [Perceptron capacity](perceptron-capacity/) | What can one neuron with plastic synapses learn? | XOR is unlearnable. Random patterns are learnable up to α = M/N ≈ 2, consistent with Cover's theorem. |
| [Finite-precision numerics](finite-precision-numerics/) | How small should a finite-difference step be? | The optimal step follows ε^(1/3) (first derivative) and ε^(1/4) (second derivative) in both single and double precision. |

<p>
<img src="neuron-models/figures/fi_curves.png" width="49%">
<img src="gene-regulation/figures/autoregulation.png" width="49%">
</p>

## Running

Python projects need only NumPy and Matplotlib:

```bash
pip install -r requirements.txt
cd neuron-models && python run.py
```

Every script uses a fixed random seed and writes its figures to `figures/`. The
numerical-analysis project is written in Wolfram Language (Mathematica).

## Background

These projects grew out of coursework at Stony Brook University: *Computational
Modeling of Physiological Systems* (AMS/BIO 332, Spring 2026) and *Foundations of
Numerical Methods* (AMS 527, Fall 2026). They were restructured into standalone
studies and extended with additional analyses and theoretical checks. Each
project's README says what was added.

See also: [nullcline-guided symbolic regression](https://github.com/seoyeonpark927/null-cline-guided-sr),
my research project on recovering equations of fast–slow systems.
