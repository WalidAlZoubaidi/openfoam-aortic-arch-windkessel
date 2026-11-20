"""Render streamlines at a given time of the aortic arch simulation (VTK, offscreen)."""
import sys
import numpy as np
import vtk
from vtk.util import numpy_support as ns

case, time, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
W, H = int(sys.argv[4]) if len(sys.argv) > 4 else 2400, int(sys.argv[5]) if len(sys.argv) > 5 else 1642

reader = vtk.vtkPOpenFOAMReader()
reader.SetFileName(f'{case}/case.foam')
reader.SetCaseType(0)  # decomposed case
reader.UpdateInformation()
reader.EnableAllPatchArrays()
reader.SetCreateCellToPoint(1)
reader.UpdateTimeStep(time) if hasattr(reader, 'UpdateTimeStep') else None
reader.GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), time)
reader.Update()
mb = reader.GetOutput()

def find(block, name):
    it = mb.NewIterator()
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        md = it.GetCurrentMetaData()
        if md.Get(vtk.vtkCompositeDataSet.NAME()) == name:
            return it.GetCurrentDataObject()
        it.GoToNextItem()
    return None

internal = find(mb, 'internalMesh')
wall = find(mb, 'wall')
inlet = find(mb, 'inlet')
print('cells', internal.GetNumberOfCells())

internal.GetPointData().SetActiveVectors('U')
umag = np.linalg.norm(ns.vtk_to_numpy(internal.GetPointData().GetArray('U')), axis=1)
print('max |U| at t =', time, ':', umag.max())

# Seeds: points spread over the inlet disc, slightly inside the domain
ib = inlet.GetBounds()
cx, cy = (ib[0] + ib[1]) / 2, (ib[2] + ib[3]) / 2
r = 0.92 * min(ib[1] - ib[0], ib[3] - ib[2]) / 2
rng = np.random.default_rng(3)
n = int(sys.argv[6]) if len(sys.argv) > 6 else 320
rr = r * np.sqrt(rng.random(n))
th = 2 * np.pi * rng.random(n)
seeds_np = np.stack([cx + rr * np.cos(th), cy + rr * np.sin(th), np.full(n, ib[4] + 0.0015)], axis=1)
seed_pts = vtk.vtkPoints()
seed_pts.SetData(ns.numpy_to_vtk(seeds_np, deep=True))
seeds = vtk.vtkPolyData()
seeds.SetPoints(seed_pts)

tracer = vtk.vtkStreamTracer()
tracer.SetInputData(internal)
tracer.SetSourceData(seeds)
tracer.SetIntegratorTypeToRungeKutta45()
tracer.SetIntegrationDirectionToForward()
tracer.SetMaximumPropagation(0.6)
tracer.SetInitialIntegrationStep(0.2)
tracer.SetMaximumNumberOfSteps(20000)
tracer.SetTerminalSpeed(1e-4)
tracer.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, 'U')
tracer.Update()

calc = vtk.vtkArrayCalculator()
calc.SetInputConnection(tracer.GetOutputPort())
calc.AddVectorArrayName('U')
calc.SetFunction('mag(U)')
calc.SetResultArrayName('Umag')
calc.Update()

# Split streamlines wherever two consecutive points are unrealistically far apart
lines_in = calc.GetOutput()
P = ns.vtk_to_numpy(lines_in.GetPoints().GetData())
conn = ns.vtk_to_numpy(lines_in.GetLines().GetData())
new_lines = vtk.vtkCellArray()
i = 0
while i < len(conn):
    npt = conn[i]
    ids = conn[i + 1:i + 1 + npt]
    i += 1 + npt
    jumps = np.where(np.linalg.norm(np.diff(P[ids], axis=0), axis=1) > 0.002)[0]
    for part in np.split(ids, jumps + 1):
        if len(part) > 20:
            new_lines.InsertNextCell(len(part))
            for pid in part:
                new_lines.InsertCellPoint(int(pid))
clean_lines = vtk.vtkPolyData()
clean_lines.SetPoints(lines_in.GetPoints())
clean_lines.SetLines(new_lines)
clean_lines.GetPointData().ShallowCopy(lines_in.GetPointData())

tubes = vtk.vtkTubeFilter()
tubes.SetInputData(clean_lines)
tubes.SetRadius(0.00013)
tubes.SetNumberOfSides(10)
tubes.CappingOn()

vmax = 1.5
lut = vtk.vtkColorTransferFunction()
# blue -> cyan -> yellow -> red (perceptually ordered, warm = fast)
for v, c in [(0.0, (0.19, 0.25, 0.62)), (0.3, (0.13, 0.56, 0.80)), (0.55, (0.36, 0.78, 0.62)),
             (0.75, (0.98, 0.80, 0.25)), (1.0, (0.86, 0.15, 0.25))]:
    lut.AddRGBPoint(v * vmax, *c)

sm = vtk.vtkPolyDataMapper()
sm.SetInputConnection(tubes.GetOutputPort())
sm.SetScalarModeToUsePointFieldData()
sm.SelectColorArray('Umag')
sm.SetLookupTable(lut)
sm.SetScalarRange(0, vmax)
sa = vtk.vtkActor()
sa.SetMapper(sm)
sa.GetProperty().SetSpecular(0.3)
sa.GetProperty().SetSpecularPower(30)

wsurf = vtk.vtkGeometryFilter()
wsurf.SetInputData(wall)
wn = vtk.vtkPolyDataNormals()
wn.SetInputConnection(wsurf.GetOutputPort())
wm = vtk.vtkPolyDataMapper()
wm.SetInputConnection(wn.GetOutputPort())
wm.ScalarVisibilityOff()
wa = vtk.vtkActor()
wa.SetMapper(wm)
wa.GetProperty().SetColor(0.93, 0.55, 0.60)
wa.GetProperty().SetOpacity(0.22)
wa.GetProperty().SetSpecular(0.4)

ren = vtk.vtkRenderer()
ren.SetBackground(1.0, 0.97, 0.975)
ren.SetBackground2(0.90, 0.95, 1.0)
ren.GradientBackgroundOn()
ren.AddActor(sa)
ren.AddActor(wa)
bar = vtk.vtkScalarBarActor()
bar.SetLookupTable(lut)
bar.SetTitle('|U| (m/s)')
bar.SetOrientationToHorizontal()
bar.SetNumberOfLabels(4)
bar.SetLabelFormat('%.1f')
bar.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
bar.GetPositionCoordinate().SetValue(0.06, 0.06)
bar.SetWidth(0.26)
bar.SetHeight(0.09)
bar.SetBarRatio(0.28)
bar.UnconstrainedFontSizeOn()
for tp in (bar.GetTitleTextProperty(), bar.GetLabelTextProperty()):
    tp.SetColor(0.2, 0.25, 0.33)
    tp.SetFontFamilyToArial()
    tp.ShadowOff()
    tp.BoldOff()
    tp.ItalicOff()
    tp.SetFontSize(int(H * 0.022))
ren.AddActor2D(bar)
ren.SetUseDepthPeeling(1)
ren.SetMaximumNumberOfPeels(8)

win = vtk.vtkRenderWindow()
win.SetOffScreenRendering(1)
win.SetAlphaBitPlanes(1)
win.SetMultiSamples(0)
win.AddRenderer(ren)
win.SetSize(W, H)

cam = ren.GetActiveCamera()
focal = np.array([0.0, -0.004, 0.016])
direction = np.array([1.0, -0.85, 0.5])
direction /= np.linalg.norm(direction)
cam.SetFocalPoint(*focal)
cam.SetPosition(*(focal + 0.27 * direction))
cam.SetViewUp(0, 0, 1)
cam.SetViewAngle(30)
ren.ResetCameraClippingRange()

win.Render()
w2i = vtk.vtkWindowToImageFilter()
w2i.SetInput(win)
w2i.Update()
png = vtk.vtkPNGWriter()
png.SetFileName(out)
png.SetInputConnection(w2i.GetOutputPort())
png.Write()
print('written', out)
