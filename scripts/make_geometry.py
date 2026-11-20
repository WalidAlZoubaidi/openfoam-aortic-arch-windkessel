"""Idealised aortic arch with three supra-aortic branches -> multi-region STL (metres)."""
import sys
import numpy as np
import vtk
from vtk.util import numpy_support as ns

out = sys.argv[1]

def polyline_tube(points, radii, n=80):
    """Resample a centreline (list of control points) with linearly varying radii."""
    pts = np.array(points, float)
    rad = np.array(radii, float)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, s[-1], n)
    P = np.stack([np.interp(t, s, pts[:, k]) for k in range(3)], axis=1)
    Rr = np.interp(t, s, rad)
    return P, Rr

def arc(theta0, theta1, n=60):
    th = np.radians(np.linspace(theta0, theta1, n))
    return [(30 * np.cos(a), -8 * np.sin(a), 30 * np.sin(a)) for a in th]

# Main aorta: ascending -> arch -> descending (mm)
main_pts = [(-30, 0, -60), (-30, 0, -10)] + arc(180, 0) + [(30, 0, -10), (30, 0, -110)]
main_rad = [12.5, 12.5] + list(np.linspace(12.5, 10.5, 60)) + [10.5, 9.5]
main = polyline_tube(main_pts, main_rad, 260)

def branch(theta, end, r0, r1):
    a = np.radians(theta)
    start = (30 * np.cos(a), -8 * np.sin(a), 30 * np.sin(a))
    mid = (start[0] + 0.3 * (end[0] - start[0]), start[1], start[2] + 18)
    return polyline_tube([start, mid, end], [r0, r0 * 0.95, r1], 60)

branches = [
    branch(128, (-30, -6, 85), 6.8, 5.8),   # brachiocephalic trunk
    branch(96, (-5, -8, 85), 4.6, 4.2),     # left common carotid
    branch(66, (22, -7, 85), 5.2, 4.6),     # left subclavian
]

# Sampling grid
h = 0.45
xs = np.arange(-48, 48 + h, h)
ys = np.arange(-24, 16 + h, h)
zs = np.arange(-103, 78 + h, h)
Z, Y, X = np.meshgrid(zs, ys, xs, indexing='ij')
Q = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

def tube_sdf(tube, chunk=400_000):
    P, Rr = tube
    A, B = P[:-1], P[1:]
    rA, rB = Rr[:-1], Rr[1:]
    AB = B - A
    L2 = (AB ** 2).sum(1)
    out = np.empty(len(Q))
    for i in range(0, len(Q), chunk):
        q = Q[i:i + chunk]
        best = np.full(len(q), np.inf)
        for j in range(len(A)):
            t = np.clip(((q - A[j]) @ AB[j]) / L2[j], 0, 1)
            d = np.linalg.norm(q - (A[j] + t[:, None] * AB[j]), axis=1) - (rA[j] + t * (rB[j] - rA[j]))
            np.minimum(best, d, out=best)
        out[i:i + chunk] = best
    return out

def smin(a, b, k):
    hh = np.clip(0.5 + 0.5 * (b - a) / k, 0, 1)
    return b * (1 - hh) + a * hh - k * hh * (1 - hh)

f = tube_sdf(main)
for br in branches:
    f = smin(f, tube_sdf(br), 3.0)

x, z = Q[:, 0], Q[:, 2]
f = np.maximum(f, np.where(x < 0, -50 - z, -100 - z))   # inlet / descending outlet
f = np.maximum(f, z - 75)                                # supra-aortic outlets

img = vtk.vtkImageData()
img.SetDimensions(len(xs), len(ys), len(zs))
img.SetOrigin(xs[0], ys[0], zs[0])
img.SetSpacing(h, h, h)
arr = ns.numpy_to_vtk(f.astype(np.float32), deep=True)
arr.SetName('f')
img.GetPointData().SetScalars(arr)

mc = vtk.vtkFlyingEdges3D()
mc.SetInputData(img)
mc.SetValue(0, 0.0)
mc.Update()

clean = vtk.vtkCleanPolyData()
clean.SetInputConnection(mc.GetOutputPort())
tri = vtk.vtkTriangleFilter()
tri.SetInputConnection(clean.GetOutputPort())
norms = vtk.vtkPolyDataNormals()
norms.SetInputConnection(tri.GetOutputPort())
norms.ComputeCellNormalsOn()
norms.ComputePointNormalsOff()
norms.ConsistencyOn()
norms.AutoOrientNormalsOn()
norms.SplittingOff()
norms.Update()
poly = norms.GetOutput()

pts = ns.vtk_to_numpy(poly.GetPoints().GetData())
cells = ns.vtk_to_numpy(poly.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
cn = ns.vtk_to_numpy(poly.GetCellData().GetNormals())
cent = pts[cells].mean(1)

region = np.full(len(cells), 'wall', dtype=object)
flat = np.abs(cn[:, 2]) > 0.95
region[flat & (np.abs(cent[:, 2] + 50) < 0.05) & (cent[:, 0] < 0)] = 'inlet'
region[flat & (np.abs(cent[:, 2] + 100) < 0.05)] = 'outlet_descending'
top = flat & (np.abs(cent[:, 2] - 75) < 0.05)
region[top & (cent[:, 0] < -15)] = 'outlet_brachiocephalic'
region[top & (cent[:, 0] >= -15) & (cent[:, 0] < 8)] = 'outlet_carotid'
region[top & (cent[:, 0] >= 8)] = 'outlet_subclavian'

with open(out, 'w') as fh:
    for name in ['wall', 'inlet', 'outlet_descending', 'outlet_brachiocephalic', 'outlet_carotid', 'outlet_subclavian']:
        idx = np.where(region == name)[0]
        print(name, len(idx))
        fh.write(f'solid {name}\n')
        for c in idx:
            n = cn[c]
            fh.write(f' facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n  outer loop\n')
            for v in pts[cells[c]] * 1e-3:
                fh.write(f'   vertex {v[0]:.8e} {v[1]:.8e} {v[2]:.8e}\n')
            fh.write('  endloop\n endfacet\n')
        fh.write(f'endsolid {name}\n')

# Areas of the openings (for the boundary conditions)
def area(name):
    idx = np.where(region == name)[0]
    p = pts[cells[idx]]
    return 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum()
for name in ['inlet', 'outlet_descending', 'outlet_brachiocephalic', 'outlet_carotid', 'outlet_subclavian']:
    print('area', name, round(area(name), 1), 'mm2')
