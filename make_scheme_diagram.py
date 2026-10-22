"""Supplementary Figure 7: kinetic diagram of the reduced KATP scheme and its correspondence with Figure 3 of the paper.
Run: python make_scheme_diagram.py  ->  figures/pub/FigS7_scheme_diagram.png"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

BLUE, GREEN, GREY, RED, PURP, ORG = '#3b6ea5', '#2e7d32', '#666666', '#c62828', '#6a1b9a', '#b26a00'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5})

fig = plt.figure(figsize=(7.6, 10.9))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(-12, 148); ax.axis('off')


def box(x, y, w, h, text, fc, ec, tc='black', fs=8.5, bold=False, lw=1.2, ls='-'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.25,rounding_size=1.2', fc=fc, ec=ec, lw=lw, ls=ls))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', color=tc, fontsize=fs, fontweight='bold' if bold else 'normal', linespacing=1.15)


def arr(p, q, color='black', both=True, lw=1.1, ls='-', style='-|>'):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='<|-|>' if both else style, mutation_scale=8, lw=lw, color=color, ls=ls, shrinkA=0, shrinkB=0))


# ------------------------------------------------------------------ panel A: ladder
ax.text(1, 146, 'A', fontsize=14, fontweight='bold')
ax.text(5, 146, 'Reduced concerted scheme: four ATP sites, closed ⇌ open of the whole channel', fontsize=9.5, fontweight='bold', va='center')
xs = [11 + 19 * j for j in range(5)]
yC, yO = 130, 110
for j, x in enumerate(xs):
    box(x - 4.2, yC - 3.5, 8.4, 7, f'C{j}', '#e3ecf7', BLUE, fs=9, bold=True)
    box(x - 4.2, yO - 3.5, 8.4, 7, f'O{j}', '#e0f1e1', GREEN, fs=9, bold=True)
    arr((x, yC - 3.7), (x, yO + 3.7), color=GREY)
    ax.text(x + 1.5, (yC + yO) / 2, ['L', 'L·r', 'L·r$^2$', 'L·r$^3$', 'L·r$^4$'][j], fontsize=8, color=GREY, va='center')
lab = ['4A/K', '3A/2K', '2A/3K', 'A/4K']
for j in range(4):
    x0, x1 = xs[j] + 4.5, xs[j + 1] - 4.5
    arr((x0, yC), (x1, yC), color=BLUE)
    arr((x0, yO), (x1, yO), color=GREEN)
    ax.text((x0 + x1) / 2, yC + 1.8, lab[j] + '$_C$', ha='center', fontsize=7, color=BLUE)
    ax.text((x0 + x1) / 2, yO + 1.8, lab[j] + '$_O$', ha='center', fontsize=7, color=GREEN)
ax.text(1, 138.5, 'closed: no current', fontsize=8, color=BLUE, fontweight='bold')
ax.text(1, 102.5, 'open: conducting', fontsize=8, color=GREEN, fontweight='bold')
ax.text(50, 97.5, 'Horizontal: ATP (A = [ATP]) binds the closed or the open conformation; each label is the ratio of the occupancies of the two states joined.\n'
        'Vertical: the gate, with C$_j$/O$_j$ = L·r$^j$, where r = K$_O$/K$_C$ > 1 (ATP binds the closed conformation more tightly).',
        fontsize=7.3, ha='center', va='center', color=GREY, linespacing=1.3)
ax.text(50, 88,
        r'$P_O=\dfrac{\sum_j O_j}{\sum_j (C_j+O_j)}=\dfrac{1}{1+L_{\rm eff}\,[(1+A/K_C)/(1+A/K_O)]^4}$,   $L_{\rm eff}=\dfrac{L}{1+D/K_D}$,  D = [MgADP]',
        ha='center', va='center', fontsize=9)

# ------------------------------------------------------------------ panel B: two-state channel
ax.text(1, 78, 'B', fontsize=14, fontweight='bold')
ax.text(5, 78, 'What the cell model uses: N identical two-state channels', fontsize=9.5, fontweight='bold', va='center')
box(6, 61, 14, 8, 'closed\n(all C$_j$)', '#e3ecf7', BLUE, bold=True)
box(40, 61, 14, 8, 'open\n(all O$_j$)', '#e0f1e1', GREEN, bold=True)
arr((20.5, 67), (39.5, 67), color='black', both=False)
arr((39.5, 63), (20.5, 63), color='black', both=False)
ax.text(30, 69, 'k$_o$ = k$_c$·P$_O$/(1−P$_O$)', ha='center', fontsize=8)
ax.text(30, 59.5, 'k$_c$ = 1/τ$_{ch}$', ha='center', fontsize=8)
ax.text(6, 55, 'P$_O$ is taken from panel A at the current ATP and\nADP (nucleotide ratio a = D/A, A + D = 5 mM).\nτ$_{ch}$ is assumed: the scheme gives P$_O$, not the speed.', fontsize=7.5, color=GREY, va='top', linespacing=1.25)
ax.text(62, 71, 'N = 10 000 channels, k of them open:', fontsize=8.5, fontweight='bold')
for dy, t in enumerate(['closing:  k → k − Binomial(k, k$_c$Δt)', 'opening:  k → k + Binomial(N − k, k$_o$Δt)', 'conductance:  g = γ·k,  γ = 70 pS', 'mean-field version: k replaced by N·P$_O$']):
    ax.text(62, 66.5 - 4.2 * dy, '• ' + t, fontsize=7.8)

# ------------------------------------------------------------------ panel C: correspondence with Figure 3
ax.text(1, 46, 'C', fontsize=14, fontweight='bold')
ax.text(5, 46, 'Correspondence with the candidate states of Figure 3 (a proposed reading)', fontsize=9.5, fontweight='bold', va='center')
cols = [4, 34, 66]
for x, h in zip(cols, ['Figure 3 node or event', 'In the reduced scheme', 'Parameter that carries it']):
    ax.text(x, 41.5, h, fontsize=8.5, fontweight='bold')
ax.plot([2, 98], [40.2, 40.2], color='black', lw=0.8)
rows = [
    ('C$_{APO}$\nclosed, no ligand', 'C0\n(closed, no ATP)', 'L\n(P$_{O,zero}$ = 1/(1+L) = 0.53)', BLUE),
    ('C$_{ATP}$  ATP-bound closed;\nATP binding and unbinding', 'C1–C4\n(ATP on the closed gate)', 'K$_C$\n(patch IC$_{50}$ 7 µM)', BLUE),
    ('O$_{NBD}$\nopen, conducting', 'O0–O4\n(the conducting set)', 'K$_O$ > K$_C$;\nP$_O$ = sum of all O$_j$', GREEN),
    ('C$_{NBD}$, C$_{PRE}$\nMgADP, NBD dimerization,\npre-open step', 'not separate states: one shift\nof the closed/open equilibrium', 'L$_{eff}$ = L/(1 + D/K$_D$)\nK$_D$ = 100 µM (assumed)', ORG),
    ('O$_{PIP2}$\nPIP$_2$-stabilized open', 'not a separate state: a change\nof the closed/open equilibrium', 'L × 4 (PIP$_2$ depletion),\nL / 4 (Q52R-like): assumed', ORG),
    ('C$_{SU}$\nsulfonylurea/glinide-bound', 'not a state: channels removed\nfrom the pool', 'available fraction 0.1\n(90 % block): assumed', ORG),
    ('Kir6.2-R201H\n(ATP site, red in Figure 3)', 'weaker ATP binding in\nboth conformations', 'K$_C$ and K$_O$ × 298/7\n(Proks et al. 2004)', RED),
]
y = 39.5
for a, b, c, col in rows:
    h = 6.3
    ax.add_patch(Rectangle((2, y - h + 0.8), 0.9, h - 1.0, fc=col, ec='none'))
    for x, t in zip(cols, (a, b, c)):
        ax.text(x, y - h / 2 + 0.3, t, fontsize=7.7, va='center', linespacing=1.15)
    y -= h
    ax.plot([2, 98], [y + 0.4, y + 0.4], color='#cccccc', lw=0.5)
ax.text(2, y - 1.0,
        'Blue and green: states of the scheme with a counterpart in Figure 3. Orange: elements of Figure 3 that the reduced scheme\n'
        'does not represent as states and reproduces only through one parameter. Figure 3 states that the relation between O$_{NBD}$\n'
        'and O$_{PIP2}$ is not established; the scheme does not distinguish them. Nothing in this correspondence is derived from structures.',
        fontsize=7.3, color=GREY, va='top', linespacing=1.3)
os.makedirs('figures/pub', exist_ok=True)
fig.savefig('figures/pub/FigS7_scheme_diagram.png', dpi=300)
print('saved')
