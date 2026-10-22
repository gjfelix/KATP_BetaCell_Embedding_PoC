"""Draws Figure 5 and the supplementary figures from results_pub.json into figures/pub/. Run: python make_pubfigs_plot.py"""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

R = json.load(open('results_pub.json'))
os.makedirs('figures/pub', exist_ok=True)
plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'axes.spines.right': False, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                     'xtick.major.size': 3, 'ytick.major.size': 3, 'legend.frameon': False})
BLUE, ORANGE, GREY, RED = '#1f5fa8', '#d9731a', '#777777', '#b3261e'


def panel(ax, lab):
    ax.text(-0.18, 1.06, lab, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')


def mean_sd(rows, key):
    vals = []
    for r in rows:
        v = [x[key] for x in r['runs'] if x[key] is not None]
        vals.append((np.mean(v) if v else np.nan, np.std(v) if v else np.nan))
    return np.array([v[0] for v in vals]), np.array([v[1] for v in vals])


# ------------------------------------------------------------------ main Figure 5
fig = plt.figure(figsize=(7.1, 6.0))
gs = GridSpec(2, 3, figure=fig, left=0.075, right=0.985, top=0.94, bottom=0.14, wspace=0.42, hspace=0.52)

# A: patch titration
S = R['scheme']
ax = fig.add_subplot(gs[0, 0]); panel(ax, 'A')
ax.semilogx(S['atp'], S['po_wt'], color=BLUE, lw=1.4, label='wild type')
ax.semilogx(S['atp'], S['po_mut'], color=ORANGE, lw=1.4, ls='--', label='R201H-like')
for ic, col in ((S['ic50_wt'], BLUE), (S['ic50_mut'], ORANGE)):
    ax.axvline(ic, color=col, lw=0.6, ls=':')
ax.text(S['ic50_wt'] * 0.8, 0.57, f"IC50 {S['ic50_wt']:.0f} µM", color=BLUE, ha='right', fontsize=6.5)
ax.text(S['ic50_mut'] * 1.25, 0.57, f"IC50 {S['ic50_mut']:.0f} µM", color=ORANGE, ha='left', fontsize=6.5)
ax.set_xlabel('[ATP] (µM), no ADP'); ax.set_ylabel('open probability'); ax.set_xlim(0.1, 1e5); ax.set_ylim(0, 0.65)
ax.text(0.22, 0.20, 'wild type', color=BLUE, fontsize=7, ha='left')
ax.text(1500, 0.30, 'R201H-like', color=ORANGE, fontsize=7, ha='left')
ax.set_title('Reduced KATP scheme\n(patch regime)', fontsize=7.5)

# B: deterministic embedding
W = R['sweep']
ax = fig.add_subplot(gs[0, 1]); panel(ax, 'B')
ra = [r['rate'] for r in W['ref']]; rA = [r['rate'] for r in W['A']]
ax.plot(W['a'], ra, 'o-', color='k', ms=2.8, lw=0.9, label='original model')
ax.plot(W['a'], rA, 's--', color=BLUE, ms=2.6, lw=0.9, label='embedded, deterministic')
ax.set_xlabel('nucleotide ratio a = ADP/ATP'); ax.set_ylabel('spike rate (Hz)')
ax.set_ylim(-0.3, 7); ax.legend(loc='upper right', fontsize=6.5)
ax.set_title('Deterministic embedding\nreproduces the original', fontsize=7.5)

# C: traces at a = 0.7
ax0 = fig.add_subplot(gs[0, 2]); ax0.axis('off'); panel(ax0, 'C')
ax0.set_title('Stochastic gating, a = 0.7\n(original: silent)', fontsize=7.5)
sub = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[0, 2], hspace=0.12)
for k, g in enumerate(('70.0', '20.0', '7.0', '2.0')):
    a_ = fig.add_subplot(sub[k]); tr = R['traces_a07'][g]
    a_.plot(tr['t'], tr['v'], lw=0.35, color=BLUE)
    a_.set_ylim(-70, 32); a_.set_yticks([])
    a_.text(0.99, 0.80, f"N·$\\mathrm{{P_O}}$ = {tr['npo']:.0f}", transform=a_.transAxes, ha='right', fontsize=6.5)
    a_.spines['left'].set_visible(False)
    if k < 3:
        a_.set_xticks([]); a_.spines['bottom'].set_visible(False)
    else:
        a_.set_xlabel('time (s)')
a_.set_xticks([20, 30, 40])

# D: rate versus N.Po at a = 0.7
rows = R['noise']['0.7']
npo = np.array([r['npo'] for r in rows]); m, s = mean_sd(rows, 'rate')
ax = fig.add_subplot(gs[1, 0]); panel(ax, 'D')
ax.axvspan(1, 10, color='#eeeeee', lw=0)
ax.errorbar(npo, m, yerr=s, fmt='o-', color=BLUE, ms=3, lw=1, capsize=1.5, elinewidth=0.6)
ax.plot(npo, [r['rate'] for r in R['det_by_gamma']['0.7']], color='k', lw=0.9, ls='--')
ax.axvline(100, color=RED, lw=0.8, ls=':')
ax.text(110, max(m) * 0.82, '1/ε²\n(ε = 10 %)', color=RED, fontsize=6.5, va='top')
ax.text(1.2, max(m) * 0.10, 'deterministic\nembedding: silent', fontsize=6.3, va='bottom')
ax.set_xscale('log'); ax.set_xlim(1, 600); ax.set_ylim(-0.25, max(m) * 1.15)
ax.set_xlabel('open channels, N·$\\mathrm{P_O}$'); ax.set_ylabel('spike rate (Hz), a = 0.7')
ax.set_title('Noise-induced spiking\nvanishes at large N·$\\mathrm{P_O}$', fontsize=7.5)

# E: CV versus N.Po at a = 0.5
rows5 = R['noise']['0.5']
npo5 = np.array([r['npo'] for r in rows5]); mc, sc = mean_sd(rows5, 'cv'); mr, _ = mean_sd(rows5, 'rate')
ax = fig.add_subplot(gs[1, 1]); panel(ax, 'E')
ax.axvspan(1, 10, color='#eeeeee', lw=0)
ax.errorbar(npo5, mc, yerr=sc, fmt='o-', color=BLUE, ms=3, lw=1, capsize=1.5, elinewidth=0.6)
ax.set_xscale('log'); ax.set_xlim(1, 600); ax.set_ylim(0, max(mc) * 1.2)
dg = np.array([r['rate'] for r in R['det_by_gamma']['0.5']]); dev = (mr - dg) / dg * 100
ax.text(0.97, 0.93, f"mean rate within {dev.max():.0f} % of the\ndeterministic embedding", transform=ax.transAxes, ha='right', va='top', fontsize=6.5)
ax.set_xlabel('open channels, N·$\\mathrm{P_O}$'); ax.set_ylabel('spike-interval variability (CV), a = 0.5')
ax.set_title('At a = 0.5 the mean rate is\nunchanged; the timing is not', fontsize=7.5)

# F: perturbations
P = R['perturb']
ax = fig.add_subplot(gs[1, 2]); panel(ax, 'F')
x = np.arange(len(P)); wdt = 0.38
det = [p['det']['rate'] for p in P]
stm = [np.mean([r['rate'] for r in p['stoch']]) for p in P]; sts = [np.std([r['rate'] for r in p['stoch']]) for p in P]
ax.bar(x - wdt / 2, det, wdt, color=BLUE, label='deterministic')
ax.bar(x + wdt / 2, stm, wdt, yerr=sts, color=ORANGE, label='stochastic (N·$\\mathrm{P_O}$ ≈ 2–3)', error_kw=dict(lw=0.6, capsize=1.2))
labs = ['wild\ntype', 'R201H\n-like', 'Q52R\n-like', 'wild\ntype', 'PIP2\ndepl.', 'sulfonyl-\nurea']
ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=6.3)
for xi, v, sm_ in zip(x, det, stm):
    if v < 0.2:
        xx = xi if sm_ < 0.2 else xi - wdt / 2
        ax.text(xx, 0.25, 'silent', rotation=90, ha='center', va='bottom', fontsize=6, color=(GREY if sm_ < 0.2 else 'white'))
ax.text(1.0, -0.27, 'a = 0.5', transform=ax.get_xaxis_transform(), ha='center', fontsize=6.5, fontweight='bold')
ax.text(4.0, -0.27, 'a = 0.7', transform=ax.get_xaxis_transform(), ha='center', fontsize=6.5, fontweight='bold')
ax.axvline(2.5, color=GREY, lw=0.5)
ax.set_ylabel('spike rate (Hz)'); ax.legend(loc='upper left', fontsize=6, bbox_to_anchor=(0.0, 1.0))
ax.set_ylim(0, max(max(det), max(np.array(stm) + np.array(sts))) * 1.45)
ax.set_title('Perturbations move the cell\nin the expected direction', fontsize=7.5)

fig.savefig('figures/pub/Fig5_PoC_human_beta_cell.png', dpi=300)
fig.savefig('figures/pub/Fig5_PoC_human_beta_cell.pdf')
plt.close(fig)

# ------------------------------------------------------------------ supplementary figures
def save(fig, name):
    fig.savefig(f'figures/pub/{name}.png', dpi=200); plt.close(fig)

# S1 implementation check
fig, ax = plt.subplots(1, 1, figsize=(3.6, 2.8))
ax.plot([r['g'] for r in R['S1']], [r['rate'] for r in R['S1']], 'o-', color='k', ms=3)
ax.axvline(0.02, color=GREY, lw=0.7, ls=':'); ax.text(0.0203, 5.5, 'annotation in the\n.ode file: 0.02', fontsize=6.5, color=GREY)
ax.set_xlabel('g_KATP (nS/pF)'); ax.set_ylabel('spike rate (Hz)'); ax.set_title('Original model (Python implementation)', fontsize=8)
fig.tight_layout(); save(fig, 'FigS1_implementation_check')

# S2 scheme variants
a_ = np.array(S['a'])
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
ax[0].plot(a_, S['target'], 'k', lw=1.2, label='original: 297 pS·a / (Nγ)')
ax[0].plot(a_, S['poA'], '--', color=BLUE, label=f"variant A (K_D = 100 µM; elasticity {S['A']['elasticity']:.2f})")
ax[0].plot(a_, S['poB'], ':', color=ORANGE, lw=2, label=f"variant B (K_D = {S['B']['KD']:.1f} µM; elasticity {S['B']['elasticity']:.2f})")
ax[0].set_xlabel('nucleotide ratio a'); ax[0].set_ylabel('open probability'); ax[0].legend(fontsize=6.2)
ax[0].set_title('Calibration of the reduced scheme', fontsize=8)
ax[1].axis('off')
txt = ('Fixed from patch data (Proks et al. 2004):\n  L = 0.89 ($\\mathrm{P_{O,zero}}$ = 0.53), patch IC50 = 7 µM\n\n'
       f"Variant A (main text):\n  K_C = {S['A']['KC']:.1f} µM, K_O = {S['A']['KO']:.0f} µM, K_D = 100 µM (assumed)\n\n"
       f"Variant B:\n  K_C = {S['B']['KC']:.1f} µM, K_O = {S['B']['KO']:.0f} µM, K_D = {S['B']['KD']:.1f} µM (fitted)\n\n"
       'Assumed: Tot = 5 mM, N = 10^4, γ = 70 pS\n(varied), τ_ch = 5 ms (varied)')
ax[1].text(0.0, 0.95, txt, va='top', fontsize=7)
fig.tight_layout(); save(fig, 'FigS2_scheme_calibration')

# S3 deterministic sweep, variants A and B
fig, ax = plt.subplots(1, 1, figsize=(4.2, 3.0))
ax.plot(W['a'], [r['rate'] for r in W['ref']], 'o-', color='k', ms=3, label='original')
ax.plot(W['a'], [r['rate'] for r in W['A']], 's--', color=BLUE, ms=3, label='embedded A')
ax.plot(W['a'], [r['rate'] for r in W['B']], '^:', color=ORANGE, ms=3, label='embedded B')
ax.set_xlabel('nucleotide ratio a'); ax.set_ylabel('spike rate (Hz)'); ax.legend(fontsize=7)
fig.tight_layout(); save(fig, 'FigS3_deterministic_sweep')

# S4 stochastic: rate, CV and silent windows versus N.Po, both a
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
for a_s, col, lab in (('0.5', BLUE, 'a = 0.5'), ('0.7', ORANGE, 'a = 0.7')):
    rows = R['noise'][a_s]; x = [r['npo'] for r in rows]
    for k, key in enumerate(('rate', 'cv', 'silent_frac')):
        m, s = mean_sd(rows, key)
        ax[k].errorbar(x, m, yerr=s, fmt='o-', color=col, ms=2.6, lw=0.9, capsize=1.2, elinewidth=0.5, label=lab)
for k, yl in enumerate(('spike rate (Hz)', 'interval CV', 'fraction of 1-s windows without spikes')):
    ax[k].set_xscale('log'); ax[k].set_xlabel('open channels, N·$\\mathrm{P_O}$'); ax[k].set_ylabel(yl, fontsize=6.8)
ax[0].legend(fontsize=6.5)
fig.tight_layout(); save(fig, 'FigS4_noise_vs_NPo')

# S5 channel kinetics
fig, ax = plt.subplots(1, 2, figsize=(6.2, 2.6))
for a_s, col in (('0.7', ORANGE), ('0.5', BLUE)):
    rows = R['tau'][a_s]; x = [r['tau'] for r in rows]
    key = 'rate' if a_s == '0.7' else 'cv'
    m, s = mean_sd(rows, key)
    ax_ = ax[0] if a_s == '0.7' else ax[1]
    ax_.errorbar(x, m, yerr=s, fmt='o-', color=col, ms=3, capsize=1.2, lw=0.9, elinewidth=0.5)
    ax_.set_xscale('log'); ax_.set_xlabel('channel relaxation time τ_ch (ms)')
ax[0].set_ylabel('spike rate (Hz), a = 0.7'); ax[1].set_ylabel('interval CV, a = 0.5')
ax[0].set_title('N·$\\mathrm{P_O}$ ≈ 3: kinetics decide the noise-induced rate', fontsize=7)
fig.tight_layout(); save(fig, 'FigS5_channel_kinetics')

# S6 perturbation traces
fig, ax = plt.subplots(len(P), 2, figsize=(7.0, 6.4), sharex=True, sharey=True)
for k, p in enumerate(P):
    for j, tr in enumerate((p['trace_det'], p['trace_stoch'])):
        ax[k, j].plot(tr['t'], tr['v'], lw=0.35, color=(BLUE, ORANGE)[j])
        ax[k, j].set_ylim(-80, 10)
    ax[k, 0].set_ylabel(f"{p['name']}\na = {p['a']}", fontsize=6.5)
ax[0, 0].set_title('deterministic embedding', fontsize=8); ax[0, 1].set_title('stochastic embedding (N·$\\mathrm{P_O}$ ≈ 2–3)', fontsize=8)
ax[-1, 0].set_xlabel('time (s)'); ax[-1, 1].set_xlabel('time (s)')
fig.tight_layout(); save(fig, 'FigS6_perturbation_traces')
print('figures written:', sorted(os.listdir('figures/pub')))
