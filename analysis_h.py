"""Spike statistics for the simulations of human_poc.py."""
import numpy as np


def spike_stats(r, t_skip=10.0, win_s=1.0):
    """Rate (Hz), ISI coefficient of variation, fraction of 1-s windows without spikes, mean V, mean open channels.
    r: result of human_poc.simulate_h. t_skip in s."""
    t = r['t']
    T = t[-1]
    sp = r['spikes_ms'] * 1e-3
    sp = sp[sp >= t_skip]
    dur = T - t_skip
    rate = len(sp) / dur
    isi = np.diff(sp)
    cv = float(isi.std() / isi.mean()) if len(isi) > 3 else float('nan')
    edges = np.arange(t_skip, T + win_s, win_s)
    cnt, _ = np.histogram(sp, bins=edges)
    silent = float((cnt == 0).mean())
    m = t >= t_skip
    return dict(rate=float(rate), cv=cv, silent_frac=silent, V_mean=float(r['V'][m].mean()),
                nopen=float(r['nopen'][m].mean()), n_spikes=int(len(sp)))


def klass(s):
    if s['rate'] < 0.2:
        return 'silent'
    if s['cv'] < 0.5:
        return 'regular spiking'
    return 'irregular/bursting'
