"""Checks the Python implementation of the original model (constant KATP conductance, no channel scheme): spike rate and
first inter-spike interval for three Heun steps and for an independent adaptive integration (LSODA)."""
import numpy as np, json, time
from scipy.integrate import solve_ivp
import human_poc as hp, analysis_h as ah

out = {}
g = 0.015
p = hp.pvec()
for dt in (0.02, 0.01, 0.005):
    r = hp.simulate_h(0, T_s=30.0, dt=dt, stride=int(round(5.0 / dt)), gk_fixed=g)
    sp = r['spikes_ms'] * 1e-3
    sp = sp[sp >= 10.0]
    out[f'Heun dt={dt} ms'] = dict(rate=len(sp) / 20.0, first_isi=float(np.diff(sp)[0]) if len(sp) > 1 else None, n=len(sp))
# independent adaptive integration (LSODA) of the same right-hand side
def f(t, y): return hp.rhs_h(y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], g, p)
y0 = [-49.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0]
t0 = time.time()
sol = solve_ivp(f, (0, 30000.0), y0, method='LSODA', rtol=1e-9, atol=1e-11, t_eval=np.arange(0, 30000.0, 0.1))
v = sol.y[0]; tt = sol.t
up = np.where((v[:-1] < -20) & (v[1:] >= -20))[0]
sp = tt[up] * 1e-3
sp = sp[sp >= 10.0]
out['LSODA rtol=1e-9'] = dict(rate=len(sp) / 20.0, first_isi=float(np.diff(sp)[0]) if len(sp) > 1 else None, n=len(sp), seconds=round(time.time() - t0, 1))
for k, v in out.items(): print(f'{k:20s}', v)
json.dump(out, open('numerics_check.json', 'w'), indent=1)
