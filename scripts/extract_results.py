"""Outlet pressures (mmHg), stroke volumes and flow split over the second heartbeat."""
import glob
import json
import numpy as np

RHO, MMHG = 1060.0, 133.322
T0, T1 = 0.8, 1.6  # second heartbeat


def load(pattern):
    data = np.loadtxt(glob.glob(pattern)[0], comments='#')
    return data[:, 0], data[:, 1]


stats = {}
for name in ['descending', 'brachiocephalic', 'carotid', 'subclavian']:
    t, p = load(f'case/postProcessing/patchAverage(patch=outlet_{name},fields=(p))/0/*.dat')
    sel = (t >= T0) & (t <= T1)
    p_mmhg = p[sel] * RHO / MMHG
    tq, q = load(f'case/postProcessing/patchFlowRate(patch=outlet_{name})/0/*.dat')
    selq = (tq >= T0) & (tq <= T1)
    volume = np.trapezoid(q[selq], tq[selq]) * 1e6
    stats[name] = {
        'systolic_mmHg': round(float(p_mmhg.max()), 1),
        'diastolic_mmHg': round(float(p_mmhg.min()), 1),
        'mean_mmHg': round(float(np.trapezoid(p_mmhg, t[sel]) / (T1 - T0)), 1),
        'volume_per_beat_mL': round(float(volume), 2),
    }
print(json.dumps(stats, indent=2))
