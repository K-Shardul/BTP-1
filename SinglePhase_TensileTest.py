import numpy as np
import os
# from sympy import Point, Line, Segment
import math
import random

os.chdir("C:\\Users\\admin\\PycharmProjects\\DualPhase")
os.getcwd()

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

L = 500.0
W = 200.0
B = 200.0
E1 = 200e09
##### GEOMETRY SPECIMEN CREATION #####
executeOnCaeStartup()
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)

stepFile = mdb.openStep('E:/IIT Bombay/UMIC/Inventor/Abaqus/DualPhase_TensileTest_200b.stp',
                    scaleFromFile=OFF)
mdb.models['Model-1'].PartFromGeometryFile(name='Part-1',
                                            geometryFile=stepFile, combine=False, dimensionality=THREE_D,
                                            type=DEFORMABLE_BODY)


# MATERIAL PROPERTY CREATION MATRIX
mdb.models['Model-1'].Material(name='First_phase')
mdb.models['Model-1'].materials['First_phase'].Elastic(table=((E1, 0.3),))


# SECTION CREATION MATRIX
mdb.models['Model-1'].HomogeneousSolidSection(name='First_phase',
    material='First_phase', thickness=None)


# SET CREATION MATRIX
p = mdb.models['Model-1'].parts['Part-1']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
p.Set(cells=cells, name='First_phase')


# SECTION ASSIGNMENT MATRIX
p = mdb.models['Model-1'].parts['Part-1']
region = p.sets['First_phase']
p = mdb.models['Model-1'].parts['Part-1']
p.SectionAssignment(region=region, sectionName='First_phase', offset=0.0,
                    offsetType=MIDDLE_SURFACE, offsetField='',
                    thicknessAssignment=FROM_SECTION)


# ASSEMBLY
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1', part=p, dependent=ON)
datum1 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=4.1405*100)
datum2 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=-4.1405*100)
p = mdb.models['Model-1'].parts['Part-1']
c = p.cells
# pickedCells = c.getSequenceFromMask(mask=('[#1 ]',), )
cells1 = c.getByBoundingBox(xMin=-4.5*100, yMin=-9.165*100, zMin=-1*100, xMax=4.5*100, yMax=9.165*100, zMax=B+1*100)
# d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=p.datums[datum1.id], cells=cells1)
p = mdb.models['Model-1'].parts['Part-1']
c = p.cells
# pickedCells = c.getSequenceFromMask(mask=('[#2 ]',), )
cells2 = c.getByBoundingBox(xMin=-4.5*100, yMin=-9.165*100, zMin=-1*100, xMax=4.5*100, yMax=6.165*100, zMax=B+1*100)
# d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=p.datums[datum2.id], cells=cells2)


# CREATING STEP
mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', timePeriod=100.0, maxNumInc=1000,
                                  initialInc=1.0, minInc=0.01, maxInc=10.0)


# I AM ABOUT TO MESH UP THE ASSEMBLY !
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
top_edges = e.getByBoundingBox(-4.5*100, 3.4255*100, -1*100, 4.5*100, 9.165*100, B+1*100)
p.seedEdgeBySize(edges=top_edges, size=50.0, deviationFactor=0.1,
                 constraint=FINER)
bottom_edges = e.getByBoundingBox(-4.5*100, -9.165*100, -1*100, 4.5*100, -3.4255*100, B+1*100)
p.seedEdgeBySize(edges=bottom_edges, size=50.0, deviationFactor=0.1,
                 constraint=FINER)
middle_edges = e.findAt(((W/2, L/2, B/4),), ((W/2, -L/4, 0),), ((W/2, -L/4, B),),
                               ((W/2, -L/2, B/4),), ((-W/2, -L/2, B/4),), ((-W/2, L/4, 0),),
                               ((-W/2, L/4, B),), ((-W/2, L/2, B/4),))
p.seedEdgeBySize(edges=middle_edges, size=10, deviationFactor=0.1,
                 constraint=FINER)
p.generateMesh()


# BOUNDARY CONDITIONS
a = mdb.models['Model-1'].rootAssembly
c1 = a.instances['Part-1'].cells
cells1 = c1.findAt(((-3.5*100, -7.165*100, B/3),))
f1 = a.instances['Part-1'].faces
faces1 = f1.findAt(((0.61874*100, -4.140554*100, B),), ((3.5*100, -7.165*100, B/3),),
                   ((1.166667*100, -8.165*100, B/3),), ((-3.5*100, -7.165*100, B/3),), ((-0.618702*100, -4.140527*100, 0.0*100),))
e1 = a.instances['Part-1'].edges
edges1 = e1.findAt(((-3.5*100, -5.915*100, B),), ((-1.75*100, -8.165*100, B),), ((3.5*100, -7.415*100, B),),
                   ((3.5*100, -7.415*100, 0.0*100),), ((3.5*100, -8.165*100, B/4),), ((-1.75*100, -8.165*100, 0.0),),
                   ((-3.5*100, -8.165*100, B/4),), ((-3.5*100, -5.915*100, 0.0),))
v1 = a.instances['Part-1'].vertices
verts1 = v1.findAt(((-3.5*100, -8.165*100, B),), ((3.5*100, -8.165*100, B),), ((3.5*100, -8.165*100, 0.0),),
                   ((-3.5*100, -8.165*100, 0.0),))
region = a.Set(vertices=verts1, edges=edges1, faces=faces1, cells=cells1,
              name='Set-1')
mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Step-1',
                                 region=region, localCsys=None)
RPid = a.ReferencePoint(point=(0.0, 8.165*100, B/2)).id
r1 = a.referencePoints
refPoints1 = (r1[RPid], )
a.Set(referencePoints=refPoints1, name='Set-Point')
region1 = a.sets['Set-Point']
mdb.models['Model-1'].DisplacementBC(name='BC-2', createStepName='Step-1',
                     region=region1, u1=0.0, u2=100.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0,
                                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='',
                                     localCsys=None)


# CONSTRAINT
s1 = a.instances['Part-1'].faces
side1Faces1 = s1.findAt(((0.0, 8.165*100, B/2),)) #, ((-3.5*100, 6.165*100, B/3),),
                        #((1.166667*100, 8.165*100, B/3),), ((3.5*100, 6.165*100, B/3),),
                        #((0.618702*100, 4.140527*100, B),))
region2 = a.Surface(side1Faces=side1Faces1, name='s_Surf-1')
mdb.models['Model-1'].Coupling(name='Constraint-1', controlPoint=region1,
                               surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                               localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)


# DEFINING OUTPUT
regionDef=mdb.models['Model-1'].rootAssembly.sets['Set-Point']
mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-2',
    createStepName='Step-1', variables=('U2', 'RF2'), region=regionDef,
    sectionPoints=DEFAULT, rebar=EXCLUDE)


# JOB SUBMISSION
a.regenerate()
mdb.Job(name='TensileTest_Spheres_Job', model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4,
        numDomains=4, numGPUs=2)
# if os.path.exists("TensileTest_Spheres_Job.lck"):
#     os.remove("TensileTest_Spheres_Job.lck")
mdb.jobs['TensileTest_Spheres_Job'].submit(consistencyChecking=OFF)


