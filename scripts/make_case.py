"""Write the OpenFOAM v14 case for the aortic arch with RCR Windkessel outlets."""
import os
import re

C = 'case'
HDR = '''FoamFile
{{
    format      ascii;
    class       {cls};
    object      {obj};
}}
'''

def write(path, obj, body, cls='dictionary'):
    os.makedirs(os.path.dirname(os.path.join(C, path)), exist_ok=True)
    with open(os.path.join(C, path), 'w') as f:
        f.write(HDR.format(cls=cls, obj=obj) + body)

# ---- split the multi-solid STL into one file per patch -------------------
stl = open(os.path.join(C, 'constant/geometry/aorta.stl')).read()
PATCHES = []
for name, content in re.findall(r'solid (\w+)\n(.*?)endsolid \1\n', stl, re.S):
    with open(os.path.join(C, f'constant/geometry/{name}.stl'), 'w') as f:
        f.write(f'solid {name}\n{content}endsolid {name}\n')
    PATCHES.append(name)
OUTLETS = [p for p in PATCHES if p.startswith('outlet')]

# ---- physiology -----------------------------------------------------------
rho = 1060.0
nu = 3.5e-3 / rho
T_beat = 0.8          # 75 bpm
T_sys = 0.3
SV = 70e-6            # stroke volume [m^3]
Q_max = SV * 3.141592653589793 / (2 * T_sys)
Q_mean = SV / T_beat
P_mean = 93 * 133.322           # [Pa]
R_tot = P_mean / Q_mean          # [Pa.s/m^3]
C_tot = 1.4e-6 / 133.322         # [m^3/Pa]
split = {
    'outlet_descending': 0.70,
    'outlet_brachiocephalic': 0.15,
    'outlet_carotid': 0.075,
    'outlet_subclavian': 0.075,
}

# ---- mesh -----------------------------------------------------------------
h = 1.0e-3
lo = (-0.050, -0.026, -0.105)
hi = (0.050, 0.018, 0.080)
n = [round((hi[i] - lo[i]) / h) for i in range(3)]
write('system/blockMeshDict', 'blockMeshDict', f'''
vertices
(
    ({lo[0]} {lo[1]} {lo[2]}) ({hi[0]} {lo[1]} {lo[2]}) ({hi[0]} {hi[1]} {lo[2]}) ({lo[0]} {hi[1]} {lo[2]})
    ({lo[0]} {lo[1]} {hi[2]}) ({hi[0]} {lo[1]} {hi[2]}) ({hi[0]} {hi[1]} {hi[2]}) ({lo[0]} {hi[1]} {hi[2]})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({n[0]} {n[1]} {n[2]}) simpleGrading (1 1 1) );
boundary ( background {{ type patch; faces ( (0 3 2 1) (4 5 6 7) (0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) ); }} );
''')

write('system/surfaceFeaturesDict', 'surfaceFeaturesDict', '''
surfaces ("aorta.stl");
#includeEtc "caseDicts/surface/surfaceFeaturesDict.cfg"
''')

geom = '\n'.join(f'    {p} {{ type triSurface; file "{p}.stl"; }}' for p in PATCHES)
refs = '\n'.join(
    f'        {p} {{ level (1 1); patchInfo {{ type {"wall" if p == "wall" else "patch"}; }} }}'
    for p in PATCHES)
write('system/snappyHexMeshDict', 'snappyHexMeshDict', f'''
#includeEtc "caseDicts/mesh/generation/snappyHexMeshDict.cfg"
castellatedMesh on;
snap            on;
addLayers       off;
geometry
{{
{geom}
}};
castellatedMeshControls
{{
    maxGlobalCells 5000000;
    nCellsBetweenLevels 3;
    features ( {{ file "aorta.eMesh"; level 1; }} );
    refinementSurfaces
    {{
{refs}
    }}
    refinementRegions {{}}
    insidePoint (-0.030 0.0 -0.030);
}}
snapControls
{{
    nSmoothPatch 5;
    tolerance 2.0;
    nSolveIter 100;
    nRelaxIter 8;
    nFeatureSnapIter 15;
    explicitFeatureSnap true;
    implicitFeatureSnap false;
}}
addLayersControls {{ layers {{}} relativeSizes true; expansionRatio 1.2; finalLayerThickness 0.5; minThickness 0.1; }}
meshQualityControls {{ }}
writeFlags ();
mergeTolerance 1e-6;
''')
write('system/meshQualityDict', 'meshQualityDict', '#includeEtc "caseDicts/mesh/generation/meshQualityDict.cfg"\n')

write('system/decomposeParDict', 'decomposeParDict', '''
numberOfSubdomains 12;
method scotch;
''')

# ---- physics --------------------------------------------------------------
write('constant/physicalProperties', 'physicalProperties', f'''
viscosityModel  constant;
nu              {nu:.6e};
''')
write('constant/momentumTransport', 'momentumTransport', 'simulationType laminar;\n')

# ---- fields ---------------------------------------------------------------
u_bc = [
    f'''    inlet
    {{
        type            flowRateInletVelocity;
        volumetricFlowRate coded;
        name            cardiacInflow;
        codeInclude
        #{{
            #include "mathematicalConstants.H"
        #}};
        code
        #{{
            // half-sine systolic ejection, {SV*1e6:.0f} mL per beat, {60/T_beat:.0f} bpm
            const scalar tc = x - {T_beat}*floor(x/{T_beat});
            return tc < {T_sys}
              ? {Q_max:.6e}*sin(constant::mathematical::pi*tc/{T_sys})
              : 0.0;
        #}};
        profile         laminarBL;
        value           uniform (0 0 0);
    }}''',
    '    wall { type noSlip; }',
] + [f'    {o} {{ type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }}' for o in OUTLETS]

write('0/U', 'U', '''
dimensions      [velocity];
internalField   uniform (0 0 0);
boundaryField
{
''' + '\n'.join(u_bc) + '\n}\n', cls='volVectorField')

def windkessel(name, frac):
    R = R_tot / frac
    Rp, Rd = 0.06 * R, 0.94 * R
    Ci = C_tot * frac
    Pc0 = Rd * Q_mean * frac
    return f'''    {name}
    {{
        type            codedFixedValue;
        value           uniform {P_mean/rho:.6f};
        name            windkessel_{name};

        code
        #{{
            // RCR parameters, SI units ({frac*100:.1f} % of cardiac output)
            const scalar Rp  = {Rp:.6e};     // [Pa.s/m^3]
            const scalar Rd  = {Rd:.6e};     // [Pa.s/m^3]
            const scalar C   = {Ci:.6e};     // [m^3/Pa]
            const scalar rho = {rho};

            static scalar Pc    = {Pc0:.6e};
            static scalar PcOld = {Pc0:.6e};
            static label  lastTimeIndex = -1;

            const label timeIndex = this->db().time().timeIndex();
            if (timeIndex != lastTimeIndex)
            {{
                PcOld = Pc;
                lastTimeIndex = timeIndex;
            }}

            const scalar Q = gSum
            (
                this->patch().lookupPatchField<surfaceScalarField, scalar>("phi")
            );

            const scalar dt = this->db().time().deltaTValue();
            Pc = (PcOld + dt/C*Q)/(1 + dt/(Rd*C));

            operator==((Rp*Q + Pc)/rho);
        #}};
    }}'''

p_bc = ['    inlet { type zeroGradient; }', '    wall { type zeroGradient; }'] + \
       [windkessel(o, split[o]) for o in OUTLETS]
write('0/p', 'p', f'''
dimensions      [kinematicPressure];
internalField   uniform {P_mean/rho:.6f};
boundaryField
{{
''' + '\n'.join(p_bc) + '\n}\n', cls='volScalarField')

# ---- numerics -------------------------------------------------------------
write('system/controlDict', 'controlDict', f'''
solver          incompressibleFluid;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {2*T_beat};
deltaT          2e-4;
writeControl    adjustableRunTime;
writeInterval   0.025;
purgeWrite      0;
writeFormat     binary;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;
adjustTimeStep  yes;
maxCo           0.9;
maxDeltaT       1e-3;

functions
{{
    #includeFunc patchFlowRate(patch=inlet)
''' + '\n'.join(f'    #includeFunc patchFlowRate(patch={o})\n    #includeFunc patchAverage(patch={o}, fields=(p))' for o in OUTLETS) + '''
}
''')

write('system/fvSchemes', 'fvSchemes', '''
ddtSchemes      { default backward; }
gradSchemes     { default Gauss linear; limited cellLimited Gauss linear 1; }
divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwind limited;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited corrected 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited corrected 0.5; }
''')

write('system/fvSolution', 'fvSolution', '''
solvers
{
    p      { solver GAMG; tolerance 1e-7; relTol 0.01; smoother DICGaussSeidel; }
    pFinal { $p; relTol 0; }
    U      { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-6; relTol 0.1; }
    UFinal { $U; relTol 0; }
}
PIMPLE
{
    nOuterCorrectors 2;
    nCorrectors      2;
    nNonOrthogonalCorrectors 1;
}
''')

print('patches', PATCHES)
print(f'Q_max={Q_max:.4e} m3/s  R_tot={R_tot:.3e}  C_tot={C_tot:.3e}  nu={nu:.3e}')
