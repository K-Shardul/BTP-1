import numpy as np
import os
# from sympy import Point, Line, Segment
import math
import random
# from scipy.optimize import minimize_scalar

os.chdir("C:\\Users\\admin\\PycharmProjects\\DualPhase")
os.getcwd()

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup


L = 500.0
W = 200.0
B = 200.0
E1 = 2e09
E2 = 200e09
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

# Distribution of Length
volume_fraction_total = 0.02
volume_fraction = 0
cylinder_count = 0
length_lognormal = []
mean = 40
sd = 8
length_upper_limit = 46
length_lower_limit = 38
r = 10

while volume_fraction < volume_fraction_total:
    fresh_lognormal_length = np.random.lognormal(mean=math.log(mean, math.exp(1)), sigma=math.log(sd, math.exp(1)), size=None)
    if fresh_lognormal_length > length_upper_limit or fresh_lognormal_length < length_lower_limit:
        continue
    volume_fraction += (math.pi * (r**2) * fresh_lognormal_length) / (L*B*W)
    cylinder_count += 1
    length_lognormal.append(fresh_lognormal_length)

    if volume_fraction > volume_fraction_total:
        volume_fraction -= (math.pi * (r**2) * fresh_lognormal_length) / (L*B*W)
        cylinder_count -= 1
        length_lognormal.pop()
        break

print(len(length_lognormal))

# Final Locations
max_incl = len(length_lognormal)
num_incl = 0
x_coordinate = []
y_coordinate = []
z_coordinate = []
long_dist = np.zeros(1000)
lat_dist_x = np.zeros(1000)
lat_dist_z = np.zeros(1000)

while num_incl < max_incl:
    random_x = random.uniform(-W / 2 + 2 * r, W / 2 - 2 * r)
    random_y = random.uniform(-L / 2 + max(length_lognormal), L / 2 - max(length_lognormal))
    random_z = random.uniform(2 * r, B - 2 * r)

    intersect_flag = False
    for j in range(0, len(x_coordinate)):
        min_long_dist = (length_lognormal[j]/2 + length_lognormal[num_incl]/2)*1.2
        min_lat_dist = 1.2*(2*r)
        long_dist[j] = abs(random_y - y_coordinate[j])
        lat_dist_x[j] = abs(random_x - x_coordinate[j])
        lat_dist_z[j] = abs(random_z - z_coordinate[j])
        lat_dist = np.sqrt((random_x - x_coordinate[j])**2 + (random_z - z_coordinate[j])**2)

        if (long_dist[j] < min_long_dist and lat_dist_x[j] < min_lat_dist) or\
                (long_dist[j] < min_long_dist and lat_dist_z[j] < min_lat_dist):
        # if lat_dist < min_lat_dist:
            intersect_flag = True
            break

    if intersect_flag == False:
        x_coordinate.append(random_x)
        y_coordinate.append(random_y)
        z_coordinate.append(random_z)
        num_incl = num_incl + 1
        print(num_incl)


# mu_A_list_all = []
# u1_list_all = []
# u2_list_all = []
# u3_list_all = []
#
#
# def P(u1, u2, u3):
#     return np.array([[u1[0], u2[0], u3[0]],
#                      [u1[1], u2[1], u3[1]],
#                      [u1[2], u2[2], u3[2]]])
#
#
# def D(w, l, b):
#     return np.array([[2 / (w ** 2), 0, 0],
#                     [0, 2 / (l ** 2), 0],
#                     [0, 0, 2 / (b ** 2)]])
#
#
# def A(P, D):
#     return np.dot(np.dot(P, D), P.T)
#
#
# def A_inv(A):
#     return np.linalg.inv(A)
#
#
# def intersection(mu_A, mu_B, Sigma_A, Sigma_B):
#     v = mu_B-mu_A
#     def k(x):
#         return 1 - np.dot(np.dot(v, np.linalg.inv((1 / (1 - x)) * Sigma_B + (1 / x) * Sigma_A)), v.T)
#     res = minimize_scalar(k, bounds=(0.001, 1), method='bounded')
#     return k(res.x) > 0
#
# while (num_incl < max_incl):
#     # random_psi = random.uniform(0, 180)
#     # random_psi_radian = random_psi * math.pi/180
#     # random_theta = random.uniform(0, 180)
#     # random_theta_radian = random_theta * math.pi / 180
#     # random_phi = random.uniform(0, 180)
#     # random_phi_radian = random_phi * math.pi / 180
#     random_x = random.uniform(-W/2+20, W/2-20)
#     random_y = random.uniform(-L/2+20, L/2-20)
#     random_z = random.uniform(20, B-20)
#     translation = [random_x, random_y, random_z]
#     l = length_lognormal[len(x_coordinate)]
#
#     mu_A = np.array(translation)
#     u1 = [1, 0, 0]
#     u2 = [0, 1, 0]
#     u3 = [0, 0, 1]
#     P_new = P(u1, u2, u3)
#     D_new = D(2*r, l, 2*r)
#     A_new = A(P_new, D_new)
#     A_inv_new = A_inv(A_new)
#
#     if num_incl == 0:
#         mu_A_list_all.append(mu_A)
#         u1_list_all.append(u1)
#         u2_list_all.append(u2)
#         u3_list_all.append(u3)
#         x_coordinate.append(random_x)
#         y_coordinate.append(random_y)
#         z_coordinate.append(random_z)
#         num_incl = num_incl + 1
#         print(num_incl)
#
#     intersect_flag = False
#     for i in range(len(mu_A_list_all)):
#         P_initial = P(u1_list_all[i], u2_list_all[i], u3_list_all[i])
#         D_initial = D(2*r, l, 2*r)
#         A_initial = A(P_initial, D_initial)
#         A_inv_initial = A_inv(A_initial)
#         if intersection(mu_A_list_all[i], mu_A, A_inv_initial, A_inv_new) == True:
#             intersect_flag = True
#             break
#
#     if intersect_flag == False:
#         mu_A_list_all.append(mu_A)
#         u1_list_all.append(u1)
#         u2_list_all.append(u2)
#         u3_list_all.append(u3)
#         x_coordinate.append(random_x)
#         y_coordinate.append(random_y)
#         z_coordinate.append(random_z)
#         num_incl = num_incl + 1
#         print(num_incl)
# while (num_incl < max_incl):
#     random_x = random.uniform(-W/2+2*length_upper_limit, W/2-2*length_upper_limit)
#     random_y = random.uniform(-L/2+2*length_upper_limit, L/2-2*length_upper_limit)
#     random_z = random.uniform(2*length_upper_limit, B-2*length_upper_limit)
#
#     isPointIntersecting = False
#     for j in range(0, len(x_coordinate)):
#         dis[j] = math.sqrt((random_x - x_coordinate[j]) ** 2 + (random_y - y_coordinate[j]) ** 2 +
#                            (random_z - z_coordinate[j]) ** 2)
#         min_dis_required = length_lognormal[num_incl] + length_lognormal[j] + 0.2
#
#         if dis[j] < min_dis_required:
#             isPointIntersecting = True
#             break
#
#     if (isPointIntersecting == False):
#         x_coordinate.append(random_x)
#         y_coordinate.append(random_y)
#         z_coordinate.append(random_z)
#         num_incl = num_incl + 1

# Material Property creation Inclusion
mdb.models['Model-1'].Material(name='Second_material')
mdb.models['Model-1'].materials['Second_material'].Elastic(table=((E2, 0.3),))

# SECTION Creation Inclusion
mdb.models['Model-1'].HomogeneousSolidSection(name='Second_material',
                                              material='Second_material', thickness=None)

# Cylinder part creation
for i in range(len(length_lognormal)):
    cylinder_name = 'Cylinder_' + str(i + 1)
    print(i)
    print(length_lognormal[i])
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
                                                sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.FixedConstraint(entity=g.findAt((0.0, 0.0)))
    s1.rectangle(point1=(0.0, length_lognormal[i]/2), point2=(r, -length_lognormal[i]/2))
    # s1.CoincidentConstraint(entity1=v.findAt((0.0, 25.0)), entity2=g.findAt((0.0,
    #                                                                          22222310.4)), addUndoState=False)
    p = mdb.models['Model-1'].Part(name=cylinder_name, dimensionality=THREE_D,
                                   type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts[cylinder_name]
    p.BaseSolidRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts[cylinder_name]

    # Set Creation Inclusion
    p1 = mdb.models['Model-1'].parts[cylinder_name]
    p = mdb.models['Model-1'].parts[cylinder_name]
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]',), )
    p.Set(cells=cells, name=cylinder_name)
    #: The set 'Sphere_1' has been created (1 cell).

    # Section Assignment Inclusion
    p = mdb.models['Model-1'].parts[cylinder_name]
    region = p.sets[cylinder_name]
    p = mdb.models['Model-1'].parts[cylinder_name]
    p.SectionAssignment(region=region, sectionName='Second_material', offset=0.0,
                        offsetType=MIDDLE_SURFACE, offsetField='',
                        thicknessAssignment=FROM_SECTION)

    # Datum Plane Partition Inclusion
    datum_yz = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0)
    datum_xy = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0)
    c = p.cells
    inclusion_cells_xy = c.findAt((0, 0.5, 0), )
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datum_xy.id], cells=inclusion_cells_xy)
    inclusion_cells_yz = c.findAt(((0, 0.5, 0.5), ), ((0, 0.5, -0.5), ))
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datum_yz.id], cells=inclusion_cells_yz)


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


# ASSEMBLY CREATION AND MOVEMENT TO DESIGNATED LOCATION
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)
for i in range(len(length_lognormal)):
    p = mdb.models['Model-1'].parts['Cylinder_'+str(i+1)]
    a.Instance(name='Cylinder_'+str(i+1)+'-1', part=p, dependent=ON)
    # print(a.Instance.)
    a = mdb.models['Model-1'].rootAssembly
    # a.rotate(instanceList=('Cylinder_' + str(i + 1) + '-1',), axisPoint=(-1.0, 0.0, 0.0),
    #          axisDirection=(2.0, 0.0, 0.0), angle=90)
    a.translate(instanceList=('Cylinder_'+str(i+1)+'-1', ), vector=(x_coordinate[i], y_coordinate[i], z_coordinate[i]))


# MERGING THE GEOMETRIES
a = mdb.models['Model-1'].rootAssembly
keys = a.instances.keys()
All_instances = []
for i in range(len(length_lognormal)+1):  # +1 coz there are len(length_lognormal) number of spheres and the Matrix part
    All_instances.append(a.instances[keys[i]])
a.InstanceFromBooleanMerge(name='Merge', instances=All_instances, keepIntersections=ON,
    originalInstances=SUPPRESS, domain=GEOMETRY)


# CREATING STEP
mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', timePeriod=100.0, maxNumInc=1000,
                                  initialInc=1.0, minInc=0.01, maxInc=10.0)


# PARTITIONING CELLS
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Merge']
a.Instance(name='Merge-1', part=p, dependent=ON)
datum1 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=4.1405*100)
datum2 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=-4.1405*100)
p = mdb.models['Model-1'].parts['Merge']
c = p.cells
# pickedCells = c.getSequenceFromMask(mask=('[#1 ]',), )
cells1 = c.getByBoundingBox(xMin=-4.5*100, yMin=-9.165*100, zMin=-1*100, xMax=4.5*100, yMax=9.165*100, zMax=B+1*100)
# d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=p.datums[datum1.id], cells=cells1)
p = mdb.models['Model-1'].parts['Merge']
c = p.cells
# pickedCells = c.getSequenceFromMask(mask=('[#2 ]',), )
cells2 = c.getByBoundingBox(xMin=-4.5*100, yMin=-9.165*100, zMin=-1*100, xMax=4.5*100, yMax=6.165*100, zMax=B+1*100)
# d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=p.datums[datum2.id], cells=cells2)


#I AM ABOUT TO MESH THINGS UP #
p = mdb.models['Model-1'].parts['Merge']
c = p.cells
e = p.edges
All_edges = e.getByBoundingBox(-1, -1, -1, 101, 101, 101)
top_edges = e.getByBoundingBox(-4.5*100, 3.4255*100, -1*100, 4.5*100, 9.165*100, B+1*100)
bottom_edges = e.getByBoundingBox(-4.5*100, -9.165*100, -1*100, 4.5*100, -3.4255*100, B+1*100)
matrix_middle_edges = e.findAt(((W/2, L/2, B/4),), ((W/2, -L/4, 0),), ((W/2, -L/4, B),),
                               ((W/2, -L/2, B/4),), ((-W/2, -L/2, B/4),), ((-W/2, L/4, 0),),
                               ((-W/2, L/4, B),), ((-W/2, L/2, B/4),))
inclusion_edges = e.getByBoundingBox(-W/2+length_upper_limit/2, -L/2+length_upper_limit/2, length_upper_limit/2,
                                     W/2-length_upper_limit/2, L/2-length_upper_limit/2, B-length_upper_limit/2)
# Seed edges
p.seedEdgeBySize(edges=inclusion_edges, size=min(length_lognormal)/2, deviationFactor=0.1, constraint=FINER)
p.seedEdgeBySize(edges=top_edges, size=min(L, B, W)/2, deviationFactor=0.1, constraint=FINER)
p.seedEdgeBySize(edges=bottom_edges, size=min(L, B, W)/2, deviationFactor=0.1, constraint=FINER)
p.seedEdgeBySize(edges=matrix_middle_edges, size=min(L, B, W)/4, deviationFactor=0.1, constraint=FINER)

# Set Mesh Controls
p.setMeshControls(regions=c.getByBoundingBox(-4.5*100, -9.165*100, -1*100, 4.5*100, 9.165*100, B+1*100),
                  elemShape=TET, technique=FREE)
elemType1 = mesh.ElemType(elemCode=C3D20R)
elemType2 = mesh.ElemType(elemCode=C3D15)
elemType3 = mesh.ElemType(elemCode=C3D10)
cells = c.getByBoundingBox(-4.5*100, -9.165*100, -1*100, 4.5*100, 9.165*100, B+1*100)
pickedRegions = (cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2,
    elemType3))

# Change element type to linear Tet
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD,
    secondOrderAccuracy=OFF, distortionControl=DEFAULT)
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2,
    elemType3))
p.generateMesh()


# BOUNDARY CONDITIONS
a = mdb.models['Model-1'].rootAssembly
c1 = a.instances['Merge-1'].cells
cells1 = c1.findAt(((-3.5*100, -7.165*100, B/3),))
f1 = a.instances['Merge-1'].faces
faces1 = f1.findAt(((0.61874*100, -4.140554*100, B),), ((3.5*100, -7.165*100, B/3),),
                   ((1.166667*100, -8.165*100, B/3),), ((-3.5*100, -7.165*100, B/3),), ((-0.618702*100, -4.140527*100, 0.0*100),))
e1 = a.instances['Merge-1'].edges
edges1 = e1.findAt(((-3.5*100, -5.915*100, B),), ((-1.75*100, -8.165*100, B),), ((3.5*100, -7.415*100, B),),
                   ((3.5*100, -7.415*100, 0.0*100),), ((3.5*100, -8.165*100, B/4),), ((-1.75*100, -8.165*100, 0.0),),
                   ((-3.5*100, -8.165*100, B/4),), ((-3.5*100, -5.915*100, 0.0),))
v1 = a.instances['Merge-1'].vertices
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
s1 = a.instances['Merge-1'].faces
side1Faces1 = s1.findAt(((0.618702*100, 4.140527*100, 0.0),), ((-3.5*100, 6.165*100, B/3),),
                        ((1.166667*100, 8.165*100, B/3),), ((3.5*100, 6.165*100, B/3),),
                        ((0.618702*100, 4.140527*100, B),))
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
mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4,
        numDomains=4, numGPUs=2)
# if os.path.exists("Job-1.lck"):
#     os.remove("Job-1.lck")
# mdb.jobs['Job-1'].submit(consistencyChecking=OFF)


# ANALYTICAL EFFECTIVE MODULUS
E_vert = 0
for i in range(len(length_lognormal)):
    E_vert += L*E1*E2/((L-length_lognormal[i])*E2 + length_lognormal[i]*E1)
E_net = (np.pi)*(r**2)*E_vert/(B*W) + (1 - (len(length_lognormal)*(np.pi)*(r**2)/(B*W)))*E1
print('E_net = ', E_net/1e09)
