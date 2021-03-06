  # This example assumes we have a mesh object selected
import os
import sys
import bpy
import csv
import math
import random
import time
import traceback
import numpy as np
import mathutils

#import matplotlib



#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d


import bmesh

#here:

#add list of all vertices of all lines
'''
================== Functions for creating the midline ===================
'''


'''
=========================================================================
SmoothFaces normal smooth
=========================================================================
'''
def SmoothFaces(obj):
    obj.data.use_auto_smooth = True
    obj.data.auto_smooth_angle = 60
    for poly in obj.data.polygons:
      poly.use_smooth = True


'''
#=============================================================================
# joinObjectsInArray
#=============================================================================
'''
def joinObjectsInArray(arr):
    if len(arr) <= 1:
        #print("joinObjectsInArray - no items")
        return False
    scene = bpy.context.scene #bpy.data.scenes['Scene']
    ctx = bpy.context.copy()
    # one of the objects to join
    ctx['active_object'] = arr[0]
    ctx['selected_objects'] = arr
    # we need the scene bases as well for joining
    ctx['selected_editable_bases'] = [scene.object_bases[ob.name] for ob in arr]
    bpy.ops.object.join(ctx)
    return arr[0]

'''
=========================================================================
DrawSphere  radius=0.1, depth=2
obejcts are created at the position of the Cursor. if it is not in 0,0,0 then
the obejcts won't be created in the way we want.
=========================================================================
'''
def DrawArrow(x,y,z, len, raduis, r,b,g, index):
    bpy.ops.mesh.primitive_cone_add(vertices=64, radius1=raduis*2, radius2=0, depth=len*0.3, view_align=False, enter_editmode=False, location=(0, 0, len*0.5 + len*0.15), layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
    cone = bpy.context.object

    bpy.ops.mesh.primitive_cylinder_add(vertices=64, radius=raduis, depth=len, view_align=False, enter_editmode=False, location=(0, 0, 0), layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
    cyll = bpy.context.object

    objToJoin = []
    objToJoin.append(cyll)
    objToJoin.append(cone)
    arrow = joinObjectsInArray(objToJoin)

    arrow.location = [0,0,0]
    for vert in arrow.data.vertices:
      vert.co.z += len/2

    arrow.rotation_euler = [x,y,z]

    mat = CreateMaterial(r,b,g)
    arrow.active_material = mat
    arrow.material_slots[0].material = mat
    arrow.material_slots[0].link = 'DATA'

    return arrow

'''
=========================================================================
DrawSphere
=========================================================================
'''
def DrawSphere():
    bpy.ops.mesh.primitive_uv_sphere_add(size=1.0, segments=64, ring_count=32, view_align=False, enter_editmode=False, location=(0, 0, 0), layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
    sphere = bpy.context.object
    SmoothFaces(sphere)

    bpy.ops.object.modifier_add(type='SUBSURF')
    sphere.modifiers["Subsurf"].show_only_control_edges = True
    sphere.modifiers["Subsurf"].levels = 1
    sphere.modifiers["Subsurf"].render_levels = 1               #increase that for quality


'''
=========================================================================
CreateMaterial
=========================================================================
'''
def CreateMaterial(r,g,b):
    mat = bpy.data.materials.new("mat")
    #mat.diffuse_intensity = 1
    mat.diffuse_color = (r,g,b,1)
    return mat
'''
=========================================================================
CreateBasis
=========================================================================
'''
def CreateBasis(points, index, thickness, curve):
    rad = thickness
    len = thickness * 20
    arZ = DrawArrow(0,0,0,len,rad, 1,0,0,0)
    arY = DrawArrow(0,1.5708,0,len,rad, 0,1,0,1)
    arX = DrawArrow(0,1.5708,1.5708,len,rad, 0,0,1,2)

    objToJoin = []
    objToJoin.append(arZ)
    objToJoin.append(arY)
    objToJoin.append(arX)
    basis = joinObjectsInArray(objToJoin)
    SmoothFaces(basis)
    basis.data.auto_smooth_angle = 1.0472			#polygon shading angle 60degrees

    #tangent constraint
    followPath = basis.constraints.new('FOLLOW_PATH')
    followPath.use_curve_follow = True
    followPath.use_curve_radius = True
    followPath.forward_axis = 'FORWARD_X'
    followPath.target = curve
    #animate follow path:
    curve.data.use_path = True
    anim = curve.data.animation_data_create()
    anim.action = bpy.data.actions.new("%sACTION" % curve.data.name)
    fcu = anim.action.fcurves.new("eval_time")
    modifier = fcu.modifiers.new("GENERATOR")
    modifier.coefficients = (0, 1)            # animation from 0 frames to 100 then looped
    # normal to sphere constraint:
    dampedPath = basis.constraints.new('DAMPED_TRACK')
    dampedPath.target = bpy.data.objects['Sphere']
    dampedPath.track_axis = 'TRACK_NEGATIVE_Y'



'''
=========================================================================
CreateCurve
=========================================================================
'''
def CreateCurve(points,thickness,colr,orderu):
  # create the Curve Datablock
  curveData = bpy.data.curves.new('myCurve', type='CURVE')
  curveData.dimensions = '3D'
  curveData.resolution_u = 3
  curveData.render_resolution_u = 3
  #thickness:
  curveData.bevel_depth = thickness                  #thickness
  curveData.bevel_resolution = 3
  curveData.fill_mode = 'FULL'
  curveData.use_path_follow = True
  curveData.use_path = True
  # map points to spline
  polyline = curveData.splines.new('NURBS')
  polyline.points.add(len(points))

  for i, coord in enumerate(points):
      x,y,z = coord
      polyline.points[i].co = (x, y, z, 1)

  curveData.splines[0].use_cyclic_u = True      #make the curve cyclic!
  curveData.splines[0].use_endpoint_u = True  # includes end points
  #curveData.splines[0].use_bezier_u = True
  curveData.splines[0].order_u = orderu

  curveOBJ = bpy.data.objects.new('myCurve', curveData)

  # deleting the mysterious end point
  #curveOBJ.select_get(True)

  scn = bpy.context.scene
  #scn.objects.link(curveOBJ)
  view_layer = bpy.context.view_layer
  '''
  bpy.context.view_layer.objects.active = curveOBJ
  curveOBJ.select_set(True, view_layer=scn.view_layers[0])
  bpy.ops.object.editmode_toggle()
  bpy.ops.curve.delete(type='VERT')
  '''


  # Link light object to the active collection of current view layer,
  # so that it'll appear in the current scene.
  view_layer.active_layer_collection.collection.objects.link(curveOBJ)# And finally select it and make it active.
  #light_object.select_set(True)
  #view_layer.objects.active = light_object
  #light_object.select_set(True)

  #===== light data here ==
  '''
  light_data = bpy.data.lights.new(name="New Light", type='POINT')

  # Create new object with our light datablock.
  light_object = bpy.data.objects.new(name="New Light", object_data=light_data)

# Link light object to the active collection of current view layer,
# so that it'll appear in the current scene.
  view_layer.active_layer_collection.collection.objects.link(light_object)

# Place light to a specified location.
  light_object.location = (5.0, 5.0, 5.0)

# And finally select it and make it active.
  light_object.select_set(True)
  view_layer.objects.active = light_object
  '''

  mat = bpy.data.materials.new("matBase")
  #mat.diffuse_intensity = 1.0
  #mat.diffuse_color = color        #selected color
  #mat.specular_hardness = 10
  #mat.specular_intensity = 0.125
  mat.diffuse_color = colr
  curveOBJ.active_material = mat
  curveOBJ.material_slots[0].link = 'OBJECT'                               #link material to object not to mesh(data)
  curveOBJ.material_slots[0].material = mat

  return curveOBJ




'''
=========================================================================
ReadDatafile
=========================================================================
    '''
def ReadDataFile(FilePath):
    global fac,N,wd,tau
    N1 = N-1
    h = 1/N
    xyz = np.loadtxt(FilePath)


    bx = np.transpose(np.zeros((N+1)))
    by = np.transpose(np.zeros((N+1)))
    bz = np.transpose(np.zeros((N+1)))

    bxp = np.transpose(np.zeros(N+1))
    byp = np.transpose(np.zeros(N+1))
    bzp = np.transpose(np.zeros(N+1))

    bx[0],by[0],bz[0] = 1,0,0
    bx[N],by[N],bz[N] = -1,0,0



    bx[1:N] = xyz[np.arange(0,3*N1-2,3)]
    by[1:N] = xyz[np.arange(1,3*N1-1,3)]
    bz[1:N] = xyz[np.arange(2,3*N1,3)]

    #bx,by,bz = xyz[0:N+1,0],xyz[0:N+1,1],xyz[0:N+1,2]

    bx,by,bz = rotate(bx,by,bz)
    '''
    fig1 = plt.figure()
    ax = fig1.add_subplot(111,projection='3d')
    ax.plot(bx,by,bz)
    plt.show()
    '''

    #========================================================
    #               b'
    #========================================================

    bxp[0]   = (bx[0]+ bx[N-1])/h
    byp[0]   = (by[0]+ by[N-1])/h
    bzp[0]   = (bz[0]+ bz[N-1])/h

    ig = np.arange(1,N+1)

    bxp[ig] = (bx[ig]-bx[ig-1])/h
    byp[ig] = (by[ig]-by[ig-1])/h
    bzp[ig] = (bz[ig]-bz[ig-1])/h


    #========================================================
    #               Constructing tangent  t= bxb'
    #========================================================
    tx = np.transpose(np.zeros((N+1)))
    ty = np.transpose(np.zeros((N+1)))
    tz = np.transpose(np.zeros((N+1)))

    #== tangent ==

    tx[0:N] = np.multiply(by[0:N],bzp[0:N]) - np.multiply(bz[0:N],byp[0:N])
    ty[0:N] = np.multiply(bz[0:N],bxp[0:N]) - np.multiply(bx[0:N],bzp[0:N])
    tz[0:N] = np.multiply(bx[0:N],byp[0:N]) - np.multiply(by[0:N],bxp[0:N])

    '''
    tx[0:N-1] = np.multiply(by[0:N-1],bz[1:N]) - np.multiply(bz[0:N-1],by[1:N])
    ty[0:N-1] = np.multiply(bz[0:N-1],bx[1:N]) - np.multiply(bx[0:N-1],bz[1:N])
    tz[0:N-1] = np.multiply(bx[0:N-1],by[1:N]) - np.multiply(by[0:N-1],bx[1:N])

    tx[N] = tx[0]#np.insert(tx,N,tx[0,0])
    ty[N] = ty[0]#np.insert(ty,N,ty[0,0])
    tz[N] = tz[0]#np.insert(tz,N,tz[0,0])
    '''
    #tau = 78.2

    tx = tx/tau
    ty = ty/tau
    tz = tz/tau
    #========================================================

    #========================================================
    #           Constructing midline
    #========================================================
    #== Midline
    rx = np.transpose(np.zeros((N+1)))
    ry = np.transpose(np.zeros((N+1)))
    rz = np.transpose(np.zeros((N+1)))

    #== Edge 1
    x1 = np.transpose(np.zeros((N+1)))
    y1 = np.transpose(np.zeros((N+1)))
    z1 = np.transpose(np.zeros((N+1)))

    #== Edge2
    x2 = np.transpose(np.zeros((N+1)))
    y2 = np.transpose(np.zeros((N+1)))
    z2 = np.transpose(np.zeros((N+1)))

    for i in range(N):
        rx[i+1] = rx[i] + h*tx[i]
        ry[i+1] = ry[i] + h*ty[i]
        rz[i+1] = rz[i] + h*tz[i]


   # Transforming the midline
    rx = rx - np.mean(rx)
    ry = ry - np.mean(ry)
    rz = rz - np.mean(rz)


    #== Half width of the strip

    l1 = sum(((rx[0:N-1]-rx[1:N])**2 + (ry[0:N-1]-ry[1:N])**2 +(rz[0:N-1]-rz[1:N])**2)**.5)
    l = l1
    #wd = l1/40
    len_b =  sum(((bx[0:N-1]-bx[1:N])**2 + (by[0:N-1]-by[1:N])**2 +(bz[0:N-1]-bz[1:N])**2)**.5)

    print(l)
    rx,ry,rz = fac*rx/l,fac*ry/l,fac*rz/l

    #== Half width of the strip

    x1 = rx-wd*bx
    y1 = ry-wd*by
    z1 = rz-wd*bz

    x2 = rx+wd*bx
    y2 = ry+wd*by
    z2 = rz+wd*bz

    # === interpolation of the midline and the edges ===
    '''
    x1,y1,z1 = fun_interpolate(x1,y1,z1)
    x2,y2,z2 = fun_interpolate(x2,y2,z2)
    rx,ry,rz = fun_interpolate(rx,ry,rz)
    '''
    #=== Array for the midline ==
    points = []
    points.append((float(rx[0]),float(ry[0]),float(rz[0])))
    for i in range(len(rx)-1):
         points.append((float(rx[i+1]),float(ry[i+1]),float(rz[i+1])))



    # Array for the edges
    data =  [ [ 0 for i in range(3) ] for j in range(4*(len(x1)-1)) ]

    for i in range(len(x1)-1):
        p1 =  np.arange(4*i,4*i+4)

        data[p1[0]][0] = x1[i]
        data[p1[0]][1] = y1[i]
        data[p1[0]][2] = z1[i]

        data[p1[1]][0] = x2[i]
        data[p1[1]][1] = y2[i]
        data[p1[1]][2] = z2[i]

        data[p1[2]][0] = x2[i+1]
        data[p1[2]][1] = y2[i+1]
        data[p1[2]][2] = z2[i+1]

        data[p1[3]][0] = x1[i+1]
        data[p1[3]][1] = y1[i+1]
        data[p1[3]][2] = z1[i+1]


    return points,data,x1,y1,z1,x2,y2,z2


'''
=========================================================================
Create face index
=========================================================================
'''
def faceindex(N):
    faces  = []

    for i in range(1,N+1):
          faces.append((int(4*i-4),int(4*i-3),int(4*i-2),int(4*i-1)))
    return faces



'''
=========================================================================
Interpolate the midline and the edges for smoothness
=========================================================================
'''
def fun_interpolate(x,y,z):
    N = len(x)-1
    N2 = 999
    x1 = np.transpose(np.zeros((N2+1)))
    y1 = np.transpose(np.zeros((N2+1)))
    z1 = np.transpose(np.zeros((N2+1)))

    #==== Linear interpolation
    s1 = np.transpose(np.linspace(0,1,N+1))
    s2 = np.transpose(np.linspace(0,1,N2+1))

    x1 = np.interp(s2,s1,x)
    y1 = np.interp(s2,s1,y)
    z1 = np.interp(s2,s1,z)

    return x1,y1,z1





def rotate(bx,by,bz):
    global th
    th = np.pi/2
    c, s = np.cos(th), np.sin(th)
    R = np.array(((c, -s), (s, c)))  # rotation matrix

    N = len(bx)-1

    #=== rotating binormal such that bz[0]=1
    for i in range(N+1):

        temp = np.dot(R,np.array((by[i],bz[i])))
        by[i]=temp[0]
        bz[i]=temp[1]

    #return bx,by,bz
    '''
    #=== rotating binormal such that tx[0]=1 or -1

    th2 = np.pi/2 - math.atan(by[1]/bx[1])
    c, s = np.cos(th2), np.sin(th2)
    R = np.array(((c,  s), (-s, c)))  # rotation matrix



    #=== rotating binormal such that bz[0]=1
    for i in range(N+1):

        temp = np.dot(R,np.array((bx[i],by[i])))
        bx[i]=temp[0]
        by[i]=temp[1]
    '''
    return bx,by,bz


'''
=========================================================================
Create rulers
=========================================================================
'''

def createRulers(x1,y1,z1,x2,y2,z2,nr):
   global thickness
   N = len(x1)-1

   dr = int(N/nr)

   ind =  np.linspace(dr,N,nr,endpoint=True)

   color = (0.925, 0.941, 0.8,1)
   th_rule= thickness/5

   for i in range(0,len(ind)):
       #=== ith ruler
       points =[]

       ind2 = int(ind[i])

       x = np.linspace(x1[ind2],x2[ind2],20)
       y = np.linspace(y1[ind2],y2[ind2],20)
       z = np.linspace(z1[ind2],z2[ind2],20)

       for j in range(0,len(x)):
           points.append([x[j],y[j],z[j]])

       orderu = 5
       curve1 = CreateCurve(points,th_rule,color,orderu)
       #==== deleting the mysterious end point of the curve
       bpy.context.view_layer.objects.active = curve1
       curve1.select_set(True)
       bpy.ops.object.editmode_toggle()
        #curveOBJ.select_set(True)
        #bpy.ops.object.editmode_toggle()
       bpy.ops.curve.delete(type='VERT')
       bpy.ops.object.editmode_toggle()


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


'''
=========================================================================
       Functions ended here. Below is the main program
=========================================================================
'''

'''
=========================================================================
 Creating mobius surface
=========================================================================
'''
global N,fac,wd,th,thickness

th = np.pi/2
fac = 5#5
#
#wd =fac/150
#wd = .01#fac/60
#wd  = fac/80  for unknotted configurations
# For schematic 1, fac = 5, wd  = fac/10
# For 5pi knot , fac = 5, wd = fac/40
# For 7pi knot , fac = 5, wd = fac/70

#=== For stable unknotted configurations
#fac = 15
wd  =  fac/20
thickness = wd/15
thickness = wd/20
# for branch 1 , n=2
'''
fac = 8, wd = fac/25, thickness = wd/20
'''
# for branch 2, n = 2
'''
fac = 15, wd = fac/45, thickness = wd/20
'''


#thickness = wd/25  # for branch 2

#thickness = wd/30  # for branch 3

nr = 6    # number of rulings



branch = 1


if (branch==1 or branch==2 or branch==3):
    N=105

if (branch==4 or branch==9):
    N=180
if(branch==5):
    N=140
if(branch== 6):
    N=90
if(branch== 7 or branch==8):
    N=150

if(branch==11 or branch==14 or branch==15):
    N=170
if(branch==12):
    N=190
if(branch==13):
    N=210

N = 240#405#190#220
N1 = N-1
h  = 1/N


#path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' + str(branch) + '/'
path = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' + str(branch) + '/'
strtau =  path + 'tau_branch_' + str(branch)+ '.txt'

tau1 = np.loadtxt(strtau)



tau2 = 78.2#10*2*np.pi   # torsion value for which I want the equilibrium configuration


p1 =  find_nearest(tau1,tau2)

#p1 = len(tau1)-1

tau  =  int(np.around(tau1[p1]*(10**10),decimals=0))

print(p1)

strb = path + 'branch_' + str(branch) + '_N' + str(N) + '_tau_' + str(tau)+ '.txt'

strb = path + 'branch_15_N300_tau_781900000000.txt'

strb = path + 'branch_15_N600_tau_782000000000_unknot1.txt'

strb = path +  'branch_15_N600_tau_782000000000_knot1.txt'
strb = path +  'branch_1_N240_tau_80936500000.txt'

tau = -8.09365
#================= Symmetric knots

#  3pi
#  1. 5pi knot
#  path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch2/'
#strb = path + '5pi_knot_symmetric.txt'
# N 240
#  2. 7pi knot
#  path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch3/'
#strb = path + '7pi_knot_symmetric.txt'
#N = 336
# 3.  9pi knot
#  path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch9/'
#strb = path +  '9pi_knot_symmetric.txt'
#N = 432


#strb = path + 'branch_1_N300_tau_850000000000.txt'
#strb = path + 'branch_3_N300_tau_714000000000.txt'
#strb = path + 'branch_4_N300_tau_807000000000.txt'
#strb = path +  'branch_6_N90_tau_334000000000.txt'# 'branch_5_N140_tau_589000000000.txt'

#==== Unknotted configurations =====
#strb = path + 'branch_6_N240_tau_279000000000.txt'
#strb = path +  'branch_6_N300_tau_600000000000.txt'
#strb = path + 'branch_1_N270_tau_342000000000.txt'

#strb = path + 'branch_8_N300_tau_470000000000.txt'
#strb = path + 'branch_4_N260_tau_410000000000.txt'

#strb = path + 'branch_12_N190_tau_850000000000.txt'
#strb = path + 'branch_13_N210_tau_850000000000.txt'
#'branch_5_N140_tau_534000000000.txt'
#strb = path + 'b_branch_2_N240_symmetric.txt'# 'branch_3_N336_tau_153797061627_symmetry.txt'#'branch_3_N168_tau_154195000000_symmetry.txt'# 'b_branch_3_N336_symmetric.txt'

#strb = path +  'branch_8_N405_tau_844150000000.txt'

print(strb)
print(fac)
#cols = (0,0.976, 0.968,1)  # rgb and facealpha for the mobius surface
cols = (0.059511,0.223228,0.708376,1)    # Eliot's suggested color
cols = (0.022406,0.110864,0.708376,1) # Darker color
points,verts,x1,y1,z1,x2,y2,z2 = ReadDataFile(strb)


N = int(len(verts)/4)
faces =  faceindex(N)
 #------------------------------------------------
# Get the active mesh
mesh = bpy.data.meshes.new("surface")

# Get a BMesh representation
bm = bmesh.new()   # create an empty BMesh
bm.from_mesh(mesh)   # fill it in from a Mesh

#add a list of vertices to the mesh
for v_co in verts:
    bm.verts.new(v_co)
#add faces ( each face is a list of 4 indices for the quad (face))
bm.verts.ensure_lookup_table()
for f_idx in faces:
     bm.faces.new([bm.verts[i] for i in f_idx])

mySurfaceObject = bpy.data.objects.new("MySurface", mesh)

# add the mesh as an object into the scene with this utility module
bpy.context.view_layer.active_layer_collection.collection.objects.link(mySurfaceObject)


mat = bpy.data.materials.new("matBase")
mat.diffuse_color = cols
mySurfaceObject.active_material = mat
mySurfaceObject.material_slots[0].link = 'OBJECT'                               #link material to object not to mesh(data)
mySurfaceObject.material_slots[0].material = mat
 #=== Removing the cube object ==







#=== Set origin to the center of mass of the object
bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS', center='BOUNDS')

# Finish up, write the bmesh back to the mesh
bm.to_mesh(mesh)
bm.free()  # free and prevent further access

# smoothening surfaces
#============ Merge vertices======
bpy.context.view_layer.objects.active = mySurfaceObject
mySurfaceObject.select_set(True)
bpy.ops.object.editmode_toggle()
bpy.ops.mesh.select_all(action='SELECT')
bpy.ops.mesh.remove_doubles()

#=========== Subdivision modifier =====
#myObject here is the one I'm applying the modifier. It must be selected/active

bpy.context.view_layer.objects.active = mySurfaceObject
mySurfaceObject.select_set(True)
bpy.ops.object.modifier_add(type='SUBSURF')                  #adds the modifier to the selected obejct
mySurfaceObject.modifiers["Subdivision"].render_levels = 6 #change factor in render
mySurfaceObject.modifiers["Subdivision"].levels = 6             #change factor in viewport only


#=========== Shade smooth ======
bpy.ops.object.editmode_toggle()
bpy.ops.object.shade_smooth()
#==============================================================================================
#----------------------------------------------------

'''
========================= Surface created =====================================
'''

'''
========================= Creating midline ===================================
'''
#points1 = []
#points = ReadDataFile(userinput['filedata'])
# Curve C10     - Trivial shape

colr = (1,0,0,1)#(0.301,0.811,0.498,1)
#points = ReadDataFile_r(strb)
orderu = 6
curve1 = CreateCurve(points, thickness, colr,orderu)            #create curve



#objs = bpy.data.objects
#objs.remove(objs["Cube"], do_unlink=True)


'''
=========================================================================
 Creating rulings
=========================================================================
'''


 #================== creating rulers
nr = 12    # number of rulings
createRulers(x1,y1,z1,x2,y2,z2,nr)












#==== Camera location ===
def update_camera(camera, focus_point=mathutils.Vector((0.0, 0.0, 1.0)), distance=4.0):
    """
    Focus the camera to a focus point and place the camera at a specific distance from that
    focus point. The camera stays in a direct line with the focus point.

    :param camera: the camera object
    :type camera: bpy.types.object
    :param focus_point: the point to focus on (default=``mathutils.Vector((0.0, 0.0, 0.0))``)
    :type focus_point: mathutils.Vector
    :param distance: the distance to keep to the focus point (default=``10.0``)
    :type distance: float
    """
    looking_direction = camera.location - focus_point
    rot_quat = looking_direction.to_track_quat('Z', 'Y')

    camera.rotation_euler = rot_quat.to_euler()
    camera.location = rot_quat@mathutils.Vector((0.0, 0.0, distance))

#update_camera(bpy.data.objects['Camera'])
