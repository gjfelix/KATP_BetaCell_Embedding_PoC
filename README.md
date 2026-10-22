# A reduced KATP channel scheme inside a human β-cell model

This repository has the code behind the illustrative simulations in our article:

> Félix-Martínez GJ and Godínez-Fernández JR. *Bridging Molecular Structures and β-Cell Models: A Structure-Informed
> Multiscale Modeling Framework for β-Cell Function and Dysfunction.* Frontiers in Chemical Biology (Hypothesis & Theory). See Figure 5 and the Supplementary Material (Section 1 and Supplementary Figures 1–3).

A word about what this is. We replace the constant KATP conductance of a published human β-cell model with a small
allosteric channel scheme, to show that the embedding works and how the cell responds when the number of open channels is
small. The channel parameters are fitted or assumed, as explained below. Nothing here is derived from protein structures or
molecular simulations, and these simulations do not test the hypothesis of the article.

## The model

**Cell model.** We use the human β-cell model of Pedersen (2010), which we re-implemented in Python. The original model file
for XPPAUT, `human_beta_model.ode`, is available from the author's page: https://www.dei.unipd.it/~pedersen/pub.html. It
is not needed to run this code, because the equations and parameters are already in `human_poc.py`. The model has no Ca²⁺
or nucleotide dynamics, and its KATP conductance is a constant, 0.015 nS/pF. It stops spiking when that conductance is
above about 0.02 nS/pF.

**Channel scheme.** The KATP channel is a concerted (Monod–Wyman–Changeux) scheme with four ATP binding sites, in which
MgADP is assumed to act by shifting the closed/open equilibrium (see Enkvetchakul and Nichols, 2003). The open probability is

$$
P_{\mathrm{O}} = \frac{1}{1 + L_{\mathrm{eff}}\left[\dfrac{1 + [\mathrm{ATP}]/K_{\mathrm{C}}}{1 + [\mathrm{ATP}]/K_{\mathrm{O}}}\right]^{4}},
\qquad
L_{\mathrm{eff}} = \frac{L}{1 + [\mathrm{MgADP}]/K_{\mathrm{D}}},
\qquad K_{\mathrm{O}} > K_{\mathrm{C}},
$$

where $L$ is the closed/open equilibrium without ligands, $K_{\mathrm{C}}$ and $K_{\mathrm{O}}$ are the ATP dissociation constants
of the closed and open conformations, and $K_{\mathrm{D}}$ describes the activation by MgADP.

In the cell, each channel is a two-state system with relaxation time $\tau_{\mathrm{ch}}$:

$$
k_{\mathrm{c}} = \frac{1}{\tau_{\mathrm{ch}}}, \qquad k_{\mathrm{o}} = k_{\mathrm{c}}\,\frac{P_{\mathrm{O}}}{1 - P_{\mathrm{O}}} .
$$

The cell has $N$ channels of conductance $\gamma$, so the mean KATP conductance is $N\,\gamma\,P_{\mathrm{O}}$. Two
versions are run:

* **Deterministic embedding:** the open fraction $P$ follows $\dfrac{dP}{dt} = k_{\mathrm{o}}(1-P) - k_{\mathrm{c}}P$.
* **Stochastic embedding:** at each time step, the number of open channels changes by two binomial draws (closing of the open
  channels with probability $k_{\mathrm{c}}\,dt$, opening of the closed channels with probability $k_{\mathrm{o}}\,dt$).

**Input.** The nucleotide ratio $a = \mathrm{ADP}/\mathrm{ATP}$ is fixed in each simulation, with
$\mathrm{ATP} = \mathrm{Tot}/(1+a)$, $\mathrm{ADP} = \mathrm{Tot}\,a/(1+a)$ and $\mathrm{Tot} = 5$ mM (assumed). To compare
with the original model, we take its conductance to be proportional to $a$ and equal to the file value at $a = 0.5$:

$$
g_{\mathrm{KATP}}(a) = 0.015\ \mathrm{nS/pF} \times \frac{a}{0.5}.
$$

**Calibration.** Two numbers come from data (Proks et al., 2004): the closed/open equilibrium without ligands, from the
open probability without ligands $P_{\mathrm{O,zero}} = 0.53$,

$$
L = \frac{1 - P_{\mathrm{O,zero}}}{P_{\mathrm{O,zero}}} = 0.89,
$$

and the ATP concentration that halves the open probability in an excised patch (ADP = 0), $\mathrm{IC}_{50} = 7$ µM.
$K_{\mathrm{C}}$ and $K_{\mathrm{O}}$ are not measured separately, so they are fitted with two conditions. First, the scheme must
give the patch $\mathrm{IC}_{50}$ above. Second, in the cell, the total conductance must follow the reference line for
nucleotide ratios between 0.40 and 0.90:

$$
N\,\gamma\,P_{\mathrm{O}}(a) = 297\ \mathrm{pS} \times a .
$$

$K_{\mathrm{D}} = 100$ µM is assumed. The fit gives

$$
K_{\mathrm{C}} = 20.0\ \mu\mathrm{M}, \qquad K_{\mathrm{O}} = 382\ \mu\mathrm{M},
$$

with a patch $\mathrm{IC}_{50}$ of 7.1 µM and $P_{\mathrm{O}}$ within 12 % of the reference line. The scheme is a bit less steep
than the original model; its elasticity at $a = 0.5$ is

$$
\frac{d\ln P_{\mathrm{O}}}{d\ln a} = 0.76 ,
$$

against 1 for the original model. We also tried fitting $K_{\mathrm{D}}$ too. It follows the line more closely but needs
$K_{\mathrm{D}} = 3.8$ µM, which is lower than usual values for MgADP, so we did not use it.

**Perturbations.** R201H-like: $K_{\mathrm{C}}$ and $K_{\mathrm{O}}$ multiplied by 298/7 (Proks et al., 2004). Q52R-like: $L$ divided
by 4. PIP₂ depletion: $L$ multiplied by 4. Sulfonylurea: 90 % of the channels unavailable. The last three are assumptions.

## What is in the repository

| File | What it does |
|---|---|
| `human_poc.py` | The Pedersen model in Python, with three modes: original, deterministic embedding and stochastic embedding |
| `katp_scheme.py` | The channel scheme, the patch $\mathrm{IC}_{50}$, the calibration and the elasticity |
| `analysis_h.py` | Spike statistics: rate, CV of the inter-spike intervals, mean potential, mean number of open channels |
| `make_pubfigs_compute.py` | Runs all the simulations of the article and writes `results_pub.json` |
| `make_pubfigs_extra.py` | Adds the deterministic spike rate for each single-channel conductance to `results_pub.json` |
| `make_pubfigs_plot.py` | Draws Figure 5 and the supplementary figures from `results_pub.json` into `figures/` |
| `make_scheme_diagram.py` | Draws the kinetic diagram of the scheme (Supplementary Figure 1) |
| `check_numerics.py`, `numerics_check.json` | Checks the Python implementation against an independent integration |
| `results_pub.json` | The results used in the article, so the figures can be redrawn without rerunning anything |
| `environment.yml` | Conda environment |

The figures are not stored in the repository. The plotting scripts create them in a `figures/` folder, with the names used
in the article: `Figure_5.png` (and `.pdf`), and `Supplementary_Figure_1.png`, `Supplementary_Figure_2.png` and
`Supplementary_Figure_3.png`.

## Installation

You need Python 3.12 (numba does not build on Python 3.14 yet). The easiest way is conda, for instance with
[Miniconda](https://docs.conda.io/en/latest/miniconda.html); the environment takes under 1 GB.

```
git clone <repository URL>
cd <repository folder>

conda env create -f environment.yml
conda activate katp-poc
```

Check that it worked:

```
python -c "import numpy, scipy, matplotlib, numba; print('ok')"
```

If you prefer pip, make a virtual environment with Python 3.12 instead:

```
python3.12 -m venv .venv
source .venv/bin/activate            # on Windows: .venv\Scripts\activate
pip install numpy==2.5.* scipy==1.18.* matplotlib==3.11.* numba==0.68.*
```

We ran everything on macOS (Apple silicon) with Python 3.12.15, numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2 and numba 0.68.0.
We have not tested other systems or the pip route.

## Running the code

Run the commands from the repository folder, with the environment active. The first run of anything that uses numba is a
little slower because it compiles the kernels.

**1. Redraw the figures from the stored results (a few seconds).**

```
python make_pubfigs_plot.py
python make_scheme_diagram.py
```

**2. Check the implementation of the cell model (about 15 s).**

```
python check_numerics.py
```

For the original model with $g_{\mathrm{KATP}} = 0.015$ nS/pF, all the integrations should give a spike rate of 4.6 Hz and a first
inter-spike interval of 0.216 s. The results are saved in `numerics_check.json`.

**3. Check the calibration of the scheme (a few seconds).**

```
python -c "import make_pubfigs_compute as C, katp_scheme as kp; c = C.cal(70.0, 'A'); print(c['KC'], c['KO'], c['ic50'], kp.elasticity(c))"
```

You should get $K_{\mathrm{C}} = 20.0$ µM, $K_{\mathrm{O}} = 382$ µM, an $\mathrm{IC}_{50}$ of 7.1 µM and an elasticity of 0.76.

**4. Rerun all the simulations (about 11 minutes on one core).** This overwrites `results_pub.json`, so keep a copy if you
want to compare.

```
python make_pubfigs_compute.py       # about 10 minutes; prints its progress
python make_pubfigs_extra.py         # about 1 minute
python make_pubfigs_plot.py          # redraws the figures
```

The stochastic runs use seeds 1, 2 and 3. On our machine the new `results_pub.json` was identical to the stored one; on
other platforms or library versions the stochastic results may differ slightly.

## Using the scheme in your own model

`katp_scheme.py` does not depend on the cell model. `po_of_a(a, L, KC, KO, KD, Tot)` gives the open probability for a
nucleotide ratio `a` and an adenine pool `Tot` (µM), and `calibrate(N, gamma, ...)` gives $K_{\mathrm{C}}$ and $K_{\mathrm{O}}$ for a
number of channels and a single-channel conductance. `simulate_h` in `human_poc.py` runs the cell model.

## Limitations

* The scheme describes equilibrium open probability only. The closing rate is $1/\tau_{\mathrm{ch}}$ with $\tau_{\mathrm{ch}}$ = 5 ms
  assumed, and the amount of noise-induced spiking depends on it (Supplementary Figure 3).
* $K_{\mathrm{C}}$ and $K_{\mathrm{O}}$ are fitted to the same conductance against which the embedded model is then compared.
* With $N = 10^4$ and $\gamma$ = 70 pS (close to the 64 pS measured in human β-cells; Ashcroft et al., 1987), there are only
  about 2–3 open channels at $a = 0.5$. The mean conductance is kept fixed when $\gamma$ is varied.

## References

* Pedersen MG (2010). A biophysical model of electrical activity in human β-cells. *Biophys J* 99:3200–3207. doi: 10.1016/j.bpj.2010.09.004
* Proks P, Antcliff JF, Lippiat J, Gloyn AL, Hattersley AT, Ashcroft FM (2004). Molecular basis of Kir6.2 mutations associated with neonatal diabetes or neonatal diabetes plus neurological features. *PNAS* 101:17539–17544. doi: 10.1073/pnas.0404756101
* Enkvetchakul D, Nichols CG (2003). Gating mechanism of KATP channels. *J Gen Physiol* 122:471–480. doi: 10.1085/jgp.200308878
* Ashcroft FM, Kakei M, Kelly RP, Sutton R (1987). ATP-sensitive K+ channels in human isolated pancreatic B-cells. *FEBS Lett* 215:9–12. doi: 10.1016/0014-5793(87)80103-5

## License and citation

MIT (see `LICENSE`). The cell model is the work of Pedersen (2010); please cite it. If you use this code, please cite the article above.
