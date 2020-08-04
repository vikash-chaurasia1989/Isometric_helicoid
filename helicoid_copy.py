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

import numpy as np
from numpy import pi as PI
from numpy import concatenate as conc

#import matplotlib



#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d


import bmesh

#here:

'''
=========================================================================
     Functions ended here. Below is the main program
=========================================================================
'''


def createplots(N,n):
    global fac, a,thickness
   #==== first creating the edges and the midline
    # initialization of the coordinates

    s = np.linspace(0,1,N+1)


    z  = s
    x1 = a*np.cos(2*PI*n*s)
    y1 = a*np.sin(2*PI*n*s)

    x2 = -a*np.cos(2*PI*n*s)
    y2 = -a*np.sin(2*PI*n*s)

    z  = fac*z



     # appending the x y z coordinates in

    points1 = []
    points2 = []
    points3 = []
    for i in range(0,N+1):
        points1.append([x1[i],y1[i],z[i]])
        points2.append([x2[i],y2[i],z[i]])
        points3.append([0*x2[i],0*y2[i],z[i]])

    colr = (1,0,0,1)#(0.301,0.811,0.498,1)
    #thickness = 2*a/20;
    #curve1 = CreateCurve(points1, thickness, colr,2)

    #==== deleting the mysterious end point of the curve
    #bpy.context.view_layer.objects.active = curve1
    #curve1.select_set(True)
    #bpy.ops.object.editmode_toggle()
    #curveOBJ.select_set(True)
    #bpy.ops.object.editmode_toggle()
    #bpy.ops.curve.delete(type='VERT')

    #curve2 = CreateCurve(points2, thickness, colr,2)
    curve1 = CreateCurve(points3, thickness, colr,2)

    led = 20
    bx1   = np.linspace(x1[0],x2[0],led)
    by1   = np.linspace(y1[0],y2[0],led)
    bz1   = z[0]*np.ones(led)

    bxN   = np.linspace(x1[N],x2[N],led)
    byN   = np.linspace(y1[N],y2[N],led)
    bzN   = z[N]*np.ones(led)

    points4 = []
    points5 = []
    for i in range(0,led):
        points4.append([bx1[i],by1[i],bz1[i]])
        points5.append([bxN[i],byN[i],bzN[i]])

    #curve4 = CreateCurve(points4, thickness, colr)
    #curve5 = CreateCurve(points5, thickness, colr)



    z1 = z
    z2 = z
    #========================================================
    #===================Curve creation commplete here ========
    #========================================================

    #================== creating rulers

    createRulers(x1,y1,z1,x2,y2,z2)




    #========================================================
    #===================Now creating surface ========
    #========================================================


    cols = (0,0.976, 0.968,1)  # rgb and facealpha for the mobius surface



    verts = helisurface(x1,y1,z1,x2,y2,z2)

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
    #----------------------------------------------------

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
    mySurfaceObject.modifiers["Subdivision"].render_levels = 2 #change factor in render
    mySurfaceObject.modifiers["Subdivision"].levels = 2             #change factor in viewport only
    mySurfaceObject.modifiers["Subdivision"].subdivision_type = 'SIMPLE'


    #=========== Shade smooth ======
    bpy.ops.object.editmode_toggle()
    bpy.ops.object.shade_smooth()




    # Drawing arrows at the ends

    points =[]
    ind2 = 0
    x = np.linspace(x1[ind2],x2[ind2],20)
    y = np.linspace(y1[ind2],y2[ind2],20)
    z = np.linspace(z1[ind2],z2[ind2],20)

    for j in range(0,len(x)):
        points.append([x[j],y[j],z[j]])
    orderu = 4

    color = (0,0,0,1)
    #curve = CreateCurve(points,thickness,color,orderu)
    #==== deleting the mysterious end point of the curve
    bpy.context.view_layer.objects.active = curve1
    curve1.select_set(True)
    bpy.ops.object.editmode_toggle()
    #curveOBJ.select_set(True)
    #bpy.ops.object.editmode_toggle()
    bpy.ops.curve.delete(type='VERT')
    bpy.ops.object.editmode_toggle()


    cols = (0,0,0,1)
    thickness = 0.001
    #CreateBasis(points, 1,thickness, curve1)
    #DrawArrow(x,y,z, .2, radius, cols, 0)
    #bpy.ops.mesh.primitive_cone_add(vertices=64, radius1=radius*2, radius2=0, depth=len*0.3,   enter_editmode=False, location=(0, 0, len*0.5 + len*0.15))
    '''
    bpy.ops.mesh.primitive_cone_add(vertices=64,radius1=raduis*2, radius2=0, depth=len1*0.3,enter_editmode=False,location=(x,y,z))

    cone = bpy.context.object
    SmoothFaces(cone)
    bpy.ops.mesh.primitive_cylinder_add(vertices=64,radius=raduis, depth=len1, enter_editmode=False, location=(x,y,z))
    cyll = bpy.context.object
    SmoothFaces(cyll)
    '''



    x = x1[N]-x2[N]
    y = y1[N]-y2[N]
    z = z1[N]-z2[N]
    # DrawArrow(x,y,z, .2, radius, cols, 0)
    '''
    ========================= Surface created =====================================
    '''

    '''
    ========================= Creating midline ===================================
    '''

 #==============================================================================
 #==================creatte plot ends here======================================
 #==============================================================================






'''
========================= Below are the subourtines needed for creating mobiud surface and the midline ==============
'''
#objs = bpy.data.objects
#objs.remove(objs["Cube"], do_unlink=True)

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

update_camera(bpy.data.objects['Camera'])



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

  for ob in bpy.context.view_layer.objects:
      if ob.type == 'MESH':
          ob.select_set(True)
          bpy.context.view_layer.objects.active = ob
      else:
          ob.select_set(False)
      bpy.ops.object.join()

  #bpy.ops.object.join(arr[0],arr[1])

  return arr[0]


 # Link light object to the active collection of current view layer,
# so that it'll appear in the current scene.
'''
=========================================================================
DrawSphere  radius=0.1, depth=2
obejcts are created at the position of the Cursor. if it is not in 0,0,0 then
the obejcts won't be created in the way we want.
=========================================================================
'''
#def DrawArrow(x,y,z, len, radius, r,b,g, index):
def DrawArrow(x,y,z, len, raduis, cols, index):
  bpy.ops.mesh.primitive_cone_add(vertices=64, radius1=raduis*2, radius2=0, depth=len*0.3,   enter_editmode=False, location=(0, 0, len*0.5 + len*0.15))
  cone = bpy.context.object

  #bpy.context.object.rotation_euler[0] = 0
  #bpy.context.object.rotation_euler[1] = 1.5708
  #bpy.context.object.rotation_euler[2] = 0
  bpy.ops.mesh.primitive_cylinder_add(vertices=64, radius=raduis, depth=len,   enter_editmode=False, location=(0, 0, 0))
  cyll = bpy.context.object

  #bpy.context.object.rotation_euler[0] = 0
  #bpy.context.object.rotation_euler[1] = 1.5708
  #bpy.context.object.rotation_euler[2] = 0
  #bpy.context.view_layer.active_layer_collection.collection.objects.link(cone)
  #bpy.context.view_layer.active_layer_collection.collection.objects.link(cyll)
  bpy.context.view_layer.objects.active = cone
  #bpy.ops.object.select_all()

  objToJoin = []
  objToJoin.append(cyll)
  objToJoin.append(cone)


  #arrow = joinObjectsInArray(objToJoin)
  #arrow = bpy.ops.object.join(objToJoin)

  arrow.location = [0,0,0]
  for vert in arrow.data.vertices:
    vert.co.z += len/2

  arrow.rotation_euler = [x,y,z]

  #mat = CreateMaterial(r,b,g)
  mat = bpy.data.materials.new("matBase")
  mat.diffuse_color = cols
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
  bpy.ops.mesh.primitive_uv_sphere_add(size=1.0, segments=64, ring_count=32,   enter_editmode=False, location=(0, 0, 0))
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
  #arZ = DrawArrow(0,0,0,len,rad, 1,0,0,0)
  #arY = DrawArrow(0,1.5708,0,len,rad, 0,1,0,1)
  #arX = DrawArrow(0,1.5708,1.5708,len,rad, 0,0,1,2)
  cols = (0,0,0,1)

  arX = DrawArrow(1,0,0,len, rad, cols, 0)

  objToJoin = []
  #objToJoin.append(arZ)
  #objToJoin.append(arY)
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

  curveData.splines[0].use_cyclic_u = False#True      #make the curve cyclic!
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
Create surface for helicoid
=========================================================================
  '''
def helisurface(x1,y1,z1,x2,y2,z2):

    N = len(x1)-1
    data =  [ [ 0 for i in range(3) ] for j in range(4*N) ]

    for i in range(N):
        p1 =  np.arange(4*i,4*i+4)

        data[p1[0]][0] = x1[i]
        data[p1[0]][1] = y1[i]
        data[p1[0]][2] = z1[i]

        data[p1[1]][0] = x2[i]
        data[p1[1]][1] = y2[i]
        data[p1[1]][2] = z2[i]

        data[p1[2]][0] = x2[i+1]
        data[p1[2]][1] = y2[i+1]
        data[p1[2]][2] = z1[i+1]

        data[p1[3]][0] = x1[i+1]
        data[p1[3]][1] = y1[i+1]
        data[p1[3]][2] = z2[i+1]


    return data


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
Create rulers
=========================================================================
'''
def createRulers(x1,y1,z1,x2,y2,z2):
   global thickness,nr
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


'''
=========================================================================
Creating mobius surface
=========================================================================
'''
global N,fac,a,th,thickness,nr


fac = 5
a  = fac/10
thickness = a/20
nr = 6 # number of rulings

N =   105
n = 8.093946633549898e+00/2/np.pi    # number of turns in the helicoid


createplots(N,n)
