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
import scipy.special as sp

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
ReadDatafile for midline
=========================================================================
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
Create face index
=========================================================================
'''
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
fac = 1
wd  = fac/50

nfold = 3*2


N     = 50*nfold
N1    = N-1

#===== values for 3fold ===

if(nfold==3*2):
    A    = 4.043808837982855e+01
    B    = 1.678602921926909e+04
    tau1 = -8.09361
elif(nfold==5*2):
    A    = 1.280530390782753e+01
    B    = 1.837823905891308e+05
    tau1 = 11.98423030321598
elif(nfold==7*2):
    A    =-8.732663887726436e+01
    B    = 7.988414727595540e+05
    tau1 = 15.37971
elif(nfold==9*2):
    A    =-2.664560362414188e+02
    B    = 2.325042020521571e+06
    tau1 = 18.47472
elif(nfold==11*2):
    A    =-5.284567289865532e+02
    B    = 5.388830397947116e+06
    tau1 = 21.35104
elif(nfold==25*2):
    #A    = 3.676534655636211e+03
    #B    = 9.711831243309331e+07
    A    = -9.727745545246296e+03
    B    =  6.598490143924065e+08
    tau1 = 78.25


    N = 1000
    N1 = N-1

al1 = 2*A + 2*np.sqrt(A**2+B)
al2 = 0
al3 = -2*A +2*np.sqrt(A**2+B)

p = np.sqrt((al3-al2)/(al3+al1))
q = np.sqrt(1-al2/al3)
r = 1/2*np.sqrt(al3+al1)
print(p)
K = sp.ellipk(p);

s = np.linspace(0,nfold*K/r,N+1);

tau = tau1*nfold*K/r;

h = 1/N*nfold*K/r;



rx = np.zeros(N+1)
ry = np.zeros(N+1)
rz = np.zeros(N+1)
tx = np.zeros(N+1)
ty = np.zeros(N+1)
tz = np.zeros(N+1)
nx = np.zeros(N+1)
ny = np.zeros(N+1)
nz = np.zeros(N+1)
bx = np.zeros(N+1)
by = np.zeros(N+1)
bz = np.zeros(N+1)

sn = np.zeros(N+1)
cn = np.zeros(N+1)
dn = np.zeros(N+1)

for i in range(0,N+1):
    sn[i],cn[i], dn[i],temp = sp.ellipj(r*s[i],p)
    #print(p**2*sn[i]**2 + dn**2)

kappa = np.sqrt(al3*(1-q**2*sn**2))

#=== correcting sign of kappa in the relevant intervals of s

#=== correcting sign of kappa in the relevant intervals of s
'''
for m in range(0,nfold-1,2):
    kappa[int((2*m+1)*N/(2*nfold)+1):int((2*m+1)*N/(2*nfold) + N/(nfold))] = -kappa[int((2*m+1)*N/(2*nfold)+1):int((2*m+1)*N/(2*nfold) + N/(nfold))]

m = nfold-1
kappa[int((2*m+1)*N/(2*nfold)+1):N+1] = -kappa[int((2*m+1)*N/(2*nfold)+1):N+1]
'''


for i in range(1,int(((nfold/2-1)/2))+1):
     #ind = int(4*(i-1)*N/nfold)+1 + temp
     ind = np.arange((4*i-3)*N/nfold,(4*i-1)*N/nfold)
     kappa[int((4*i-3)*N/nfold):int((4*i-1)*N/nfold)] = -kappa[int((4*i-3)*N/nfold):int((4*i-1)*N/nfold)]

#ind = np.arange(int((nfold-1)*N/nfold),(N+1))
kappa[int((nfold-1)*N/nfold):(N+1)] = -kappa[int((nfold-1)*N/nfold):(N+1)]



# === frenet frame integration
initial = np.array([0, -0.1276, 0, 0.9280, 0, -0.3726, 0.3726, 0, 0.9280, 0, 1, 0],'d')
initial = np.array([0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0],'d')

rx[0] = initial[0]
ry[0] = initial[1]
rz[0] = initial[2]

tx[0] = initial[3]
ty[0] = initial[4]
tz[0] = initial[5]

nx[0] = initial[6]
ny[0] = initial[7]
nz[0] = initial[8]

bx[0] = initial[9]
by[0] = initial[10]
bz[0] = initial[11]

i = 0

rx[i+1] = rx[i] + h*tx[i];
ry[i+1] = ry[i] + h*ty[i];
rz[i+1] = rz[i] + h*tz[i];


tx[i+1] = tx[i] + kappa[i]*h*nx[i]
ty[i+1] = ty[i] + kappa[i]*h*ny[i]
tz[i+1] = tz[i] + kappa[i]*h*nz[i]

nx[i+1] = nx[i] + (-kappa[i]*tx[i] + tau*bx[i])*h
ny[i+1] = ny[i] + (-kappa[i]*ty[i] + tau*by[i])*h
nz[i+1] = nz[i] + (-kappa[i]*tz[i] + tau*bz[i])*h

bx[i+1] = bx[i] - tau*nx[i]*h
by[i+1] = by[i] - tau*ny[i]*h
bz[i+1] = bz[i] - tau*nz[i]*h

#== central difference integration for t,n, and b

for i in range(1,N):

    rx[i] = rx[i-1]+h*tx[i-1]
    ry[i] = ry[i-1]+h*ty[i-1]
    rz[i] = rz[i-1]+h*tz[i-1]

    tx[i+1] = tx[i-1] + kappa[i]*2*h*nx[i]
    ty[i+1] = ty[i-1] + kappa[i]*2*h*ny[i]
    tz[i+1] = tz[i-1] + kappa[i]*2*h*nz[i]

    nx[i+1] = nx[i-1] + (-kappa[i]*tx[i] + tau*bx[i])*2*h
    ny[i+1] = ny[i-1] + (-kappa[i]*ty[i] + tau*by[i])*2*h
    nz[i+1] = nz[i-1] + (-kappa[i]*tz[i] + tau*bz[i])*2*h

    bx[i+1] = bx[i-1] - tau*nx[i]*2*h
    by[i+1] = by[i-1] - tau*ny[i]*2*h
    bz[i+1] = bz[i-1] - tau*nz[i]*2*h

i = N
rx[i] = rx[i-1]+h*tx[i-1]
ry[i] = ry[i-1]+h*ty[i-1]
rz[i] = rz[i-1]+h*tz[i-1]

#===== Integration complete ===

#==== Constructing midline and edges
# Transforming the midline
rx = rx - np.mean(rx)
ry = ry - np.mean(ry)
rz = rz - np.mean(rz)


#== Half width of the strip

l = sum(((rx[0:N-1]-rx[1:N])**2 + (ry[0:N-1]-ry[1:N])**2 +(rz[0:N-1]-rz[1:N])**2)**.5)

rx,ry,rz = fac*rx/l,fac*ry/l,fac*rz/l

points = []
points.append((float(rx[0]),float(ry[0]),float(rz[0])))
for i in range(N):
     points.append((float(rx[i+1]),float(ry[i+1]),float(rz[i+1])))

#== Half width of the strip
x1 = rx-wd*bx
y1 = ry-wd*by
z1 = rz-wd*bz

x2 = rx+wd*bx
y2 = ry+wd*by
z2 = rz+wd*bz


verts =  [ [ 0 for i in range(3) ] for j in range(4*N) ]

for i in range(N):
    p1 =  np.arange(4*i,4*i+4)

    verts[p1[0]][0] = x1[i]
    verts[p1[0]][1] = y1[i]
    verts[p1[0]][2] = z1[i]

    verts[p1[1]][0] = x2[i]
    verts[p1[1]][1] = y2[i]
    verts[p1[1]][2] = z2[i]

    verts[p1[2]][0] = x2[i+1]
    verts[p1[2]][1] = y2[i+1]
    verts[p1[2]][2] = z2[i+1]

    verts[p1[3]][0] = x1[i+1]
    verts[p1[3]][1] = y1[i+1]
    verts[p1[3]][2] = z1[i+1]



thickness = wd/20


#cols = (0,0.976, 0.968,1)  # rgb and facealpha for the mobius surface
cols = (0.059511,0.223228,0.708376,1)    # Eliot's suggested color
cols = (0.022406,0.110864,0.708376,1) # Darker color

#verts,x1,y1,z1,x2,y2,z2 = ReadDataFile(strb)

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
nr = 12
#createRulers(x1,y1,z1,x2,y2,z2,nr)












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
