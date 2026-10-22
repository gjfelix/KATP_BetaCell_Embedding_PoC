"""Adds to results_pub.json the deterministic spike rate for each single-channel conductance, with the same calibration
as the stochastic runs. Run after make_pubfigs_compute.py."""
import json
import make_pubfigs_compute as C
R = json.load(open('results_pub.json'))
GAM = [140.0, 100.0, 70.0, 50.0, 35.0, 25.0, 20.0, 14.0, 10.0, 7.0, 5.0, 3.5, 2.5, 2.0, 1.0, 0.5]
R['det_by_gamma'] = {str(a): [dict(gamma=g, **C.st(C.sim(1, a=a, gamma=g, T_s=60.0))) for g in GAM] for a in (0.5, 0.7)}
json.dump(R, open('results_pub.json', 'w'))
for a in ('0.5', '0.7'):
    print(a, [(r['gamma'], round(r['rate'], 2)) for r in R['det_by_gamma'][a]])
