"""Mid-plane slice of the aortic arch: velocity magnitude at a given time, or the mesh."""
import sys
import numpy as np
import vtk

case, time, mode, out = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
W, H = 2000, 1300

reader = vtk.vtkPOpenFOAMReader()
reader.SetFileName(f'{case}/case.foam')
reader.SetCaseType(0)
reader.UpdateInformation()
reader.EnableAllPatchArrays()
reader.SetCreateCellToPoint(1)
reader.GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), time)
reader.Update()
mb = reader.GetOutput()

def find(name):
    it = mb.NewIterator()
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == name:
            return it.GetCurrentDataObject()
        it.GoToNextItem()

internal, wall = find('internalMesh'), find('wall')

plane = vtk.vtkPlane()
plane.SetOrigin(0, -0.006, 0)
plane.SetNormal(0, 1, 0)

ren = vtk.vtkRenderer()
ren.SetBackground(1, 1, 1)

if mode == 'velocity':
    cut = vtk.vtkCutter()
    cut.SetCutFunction(plane)
    cut.SetInputData(internal)
    calc = vtk.vtkArrayCalculator()
    calc.SetInputConnection(cut.GetOutputPort())
    calc.AddVectorArrayName('U')
    calc.SetFunction('mag(U)')
    calc.SetResultArrayName('Umag')
    vmax = 1.5
    lut = vtk.vtkColorTransferFunction()
    for v, c in [(0.0, (0.19, 0.25, 0.62)), (0.3, (0.13, 0.56, 0.80)), (0.55, (0.36, 0.78, 0.62)),
                 (0.75, (0.98, 0.80, 0.25)), (1.0, (0.86, 0.15, 0.25))]:
        lut.AddRGBPoint(v * vmax, *c)
    m = vtk.vtkPolyDataMapper()
    m.SetInputConnection(calc.GetOutputPort())
    m.SetScalarModeToUsePointFieldData()
    m.SelectColorArray('Umag')
    m.SetLookupTable(lut)
    m.SetScalarRange(0, vmax)
    a = vtk.vtkActor()
    a.SetMapper(m)
    a.GetProperty().LightingOff()
    ren.AddActor(a)
    bar = vtk.vtkScalarBarActor()
    bar.SetLookupTable(lut)
    bar.SetTitle('U (m/s)')
    bar.SetOrientationToHorizontal()
    bar.SetNumberOfLabels(4)
    bar.SetLabelFormat('%.1f')
    bar.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    bar.GetPositionCoordinate().SetValue(0.06, 0.05)
    bar.SetWidth(0.28)
    bar.SetHeight(0.09)
    bar.SetBarRatio(0.28)
    bar.UnconstrainedFontSizeOn()
    for tp in (bar.GetTitleTextProperty(), bar.GetLabelTextProperty()):
        tp.SetColor(0.2, 0.25, 0.33)
        tp.SetFontFamilyToArial()
        tp.ShadowOff()
        tp.BoldOff()
        tp.ItalicOff()
        tp.SetFontSize(26)
    ren.AddActor2D(bar)
    focal, dist = (0.0, -0.006, 0.0), 0.30
else:
    # Cells crossed by the plane, drawn with their edges
    ext = vtk.vtkExtractGeometry()
    ext.SetInputData(internal)
    ext.SetImplicitFunction(plane)
    ext.ExtractInsideOn()
    ext.ExtractBoundaryCellsOn()
    ext.ExtractOnlyBoundaryCellsOn()
    geo = vtk.vtkGeometryFilter()
    geo.SetInputConnection(ext.GetOutputPort())
    m = vtk.vtkPolyDataMapper()
    m.SetInputConnection(geo.GetOutputPort())
    m.ScalarVisibilityOff()
    a = vtk.vtkActor()
    a.SetMapper(m)
    p = a.GetProperty()
    p.SetColor(0.93, 0.95, 0.98)
    p.EdgeVisibilityOn()
    p.SetEdgeColor(0.25, 0.32, 0.45)
    p.SetLineWidth(1.0)
    ren.AddActor(a)
    focal, dist = (-0.012, -0.006, 0.030), 0.075

# Faint wall outline for context
ws = vtk.vtkGeometryFilter()
ws.SetInputData(wall)
wm = vtk.vtkPolyDataMapper()
wm.SetInputConnection(ws.GetOutputPort())
wm.ScalarVisibilityOff()
wa = vtk.vtkActor()
wa.SetMapper(wm)
wa.GetProperty().SetColor(0.85, 0.55, 0.6)
wa.GetProperty().SetOpacity(0.08)
ren.AddActor(wa)

win = vtk.vtkRenderWindow()
win.SetOffScreenRendering(1)
win.AddRenderer(ren)
win.SetSize(W, H)
cam = ren.GetActiveCamera()
cam.SetFocalPoint(*focal)
cam.SetPosition(focal[0], focal[1] - dist, focal[2])
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
