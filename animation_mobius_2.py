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
ReadDatafile for midline
=========================================================================
'''
def ReadDataFile_r(FilePath):
    points = []
    #with open(FilePath,"r",encoding='utf-8', errors='ignore') as f:
    with open(FilePath,'r') as f:

     # c = csv.reader(f, delimiter=',', skipinitialspace=True)
      for sRow in f:
          a,b,c = sRow.split()                                                           # thats an individual object
          points.append((float(a),float(b),float(c)))

      N = int(len(points)-1)
      xyz = np.array(points)

      bx = np.transpose(np.zeros((N+1)))
      by = np.transpose(np.zeros((N+1)))
      bz = np.transpose(np.zeros((N+1)))

      bx,by,bz = xyz[0:N+1,0],xyz[0:N+1,1],xyz[0:N+1,2]

      bxp = np.transpose(np.zeros(N+1))
      byp = np.transpose(np.zeros(N+1))
      bzp = np.transpose(np.zeros(N+1))

      bx2p = np.transpose(np.zeros(N+1))
      by2p = np.transpose(np.zeros(N+1))
      bz2p = np.transpose(np.zeros(N+1))


      #========================================================
      #               b'
      #========================================================
      h = 1/N
      bxp[0]   = (bx[0]+ bx[N-1])/h
      byp[0]   = (by[0]+ by[N-1])/h
      bzp[0]   = (bz[0]+ bz[N-1])/h

      bx2p[0]  = (bx[1] - bx[N-1] -2*bx[0])/h**2
      by2p[0]  = (by[1] - by[N-1] -2*by[0])/h**2
      bz2p[0]  = (bz[1] - bz[N-1] -2*bz[0])/h**2



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

      #========================================================

      #========================================================
      #           Constructing midline
      #========================================================
      #== Midline
      rx = np.transpose(np.zeros((N+1)))
      ry = np.transpose(np.zeros((N+1)))
      rz = np.transpose(np.zeros((N+1)))

      #== Edge 1
      rx1 = np.transpose(np.zeros((N+1)))
      ry1 = np.transpose(np.zeros((N+1)))
      rz1 = np.transpose(np.zeros((N+1)))

      #== Edge2
      rx2 = np.transpose(np.zeros((N+1)))
      ry2 = np.transpose(np.zeros((N+1)))
      rz2 = np.transpose(np.zeros((N+1)))

      for i in range(N):
          rx[i+1] = rx[i] + tx[i]
          ry[i+1] = ry[i] + ty[i]
          rz[i+1] = rz[i] + tz[i]


     # Transforming the midline
      rx = rx - np.mean(rx)
      ry = ry - np.mean(ry)
      rz = rz - np.mean(rz)


      #== Half width of the strip

      l = sum(((rx[0:N-1]-rx[1:N])**2 + (ry[0:N-1]-ry[1:N])**2 +(rz[0:N-1]-rz[1:N])**2)**.5)

      rx,ry,rz = rx/l,ry/l,rz/l

      points = []
      points.append((float(rx[0]),float(ry[0]),float(rz[0])))
      for i in range(N):
           points.append((float(rx[i+1]),float(ry[i+1]),float(rz[i+1])))
    return points



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
def CreateCurve(points,thickness,colr):
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

    curveOBJ = bpy.data.objects.new('myCurve', curveData)
    scn = bpy.context.scene
    #scn.objects.link(curveOBJ)
    view_layer = bpy.context.view_layer
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
ReadDatafile
=========================================================================
    '''
def ReadDataFile(FilePath):
    global wd
    points = []
    with open(FilePath,'r') as f:
#      c = csv.reader(f, delimiter=',', skipinitialspace=True)
      # c = csv.reader(f, delimiter=',', skipinitialspace=True)
      for sRow in f:
          a,b,c = sRow.split()                                                                 # thats an individual object
          points.append((float(a),float(b),float(c)))
      #==== Converting from object to numerical array

      N = int(len(points)-1)
      xyz = np.array(points)

      bx = np.transpose(np.zeros((N+1)))
      by = np.transpose(np.zeros((N+1)))
      bz = np.transpose(np.zeros((N+1)))

      bx,by,bz = xyz[0:N+1,0],xyz[0:N+1,1],xyz[0:N+1,2]

      bxp = np.transpose(np.zeros(N+1))
      byp = np.transpose(np.zeros(N+1))
      bzp = np.transpose(np.zeros(N+1))

      bx2p = np.transpose(np.zeros(N+1))
      by2p = np.transpose(np.zeros(N+1))
      bz2p = np.transpose(np.zeros(N+1))


      #========================================================
      #               b'
      #========================================================
      h = 1/N
      bxp[0]   = (bx[0]+ bx[N-1])/h
      byp[0]   = (by[0]+ by[N-1])/h
      bzp[0]   = (bz[0]+ bz[N-1])/h

      bx2p[0]  = (bx[1] - bx[N-1] -2*bx[0])/h**2
      by2p[0]  = (by[1] - by[N-1] -2*by[0])/h**2
      bz2p[0]  = (bz[1] - bz[N-1] -2*bz[0])/h**2



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

      #========================================================

      #========================================================
      #           Constructing midline
      #========================================================
      #== Midline
      rx = np.transpose(np.zeros((N+1)))
      ry = np.transpose(np.zeros((N+1)))
      rz = np.transpose(np.zeros((N+1)))

      #== Edge 1
      rx1 = np.transpose(np.zeros((N+1)))
      ry1 = np.transpose(np.zeros((N+1)))
      rz1 = np.transpose(np.zeros((N+1)))

      #== Edge2
      rx2 = np.transpose(np.zeros((N+1)))
      ry2 = np.transpose(np.zeros((N+1)))
      rz2 = np.transpose(np.zeros((N+1)))

      for i in range(N):
          rx[i+1] = rx[i] + tx[i]
          ry[i+1] = ry[i] + ty[i]
          rz[i+1] = rz[i] + tz[i]


     # Transforming the midline
      rx = rx - np.mean(rx)
      ry = ry - np.mean(ry)
      rz = rz - np.mean(rz)


      #== Half width of the strip

      l = sum(((rx[0:N-1]-rx[1:N])**2 + (ry[0:N-1]-ry[1:N])**2 +(rz[0:N-1]-rz[1:N])**2)**.5)

      rx,ry,rz = rx/l,ry/l,rz/l

      #==== Transforming the midline such that its symmetry plane is parallel to z axis



      #== Half width of the strip




      rx1 = rx-wd*bx
      ry1 = ry-wd*by
      rz1 = rz-wd*bz

      rx2 = rx+wd*bx
      ry2 = ry+wd*by
      rz2 = rz+wd*bz


      data =  [ [ 0 for i in range(3) ] for j in range(4*N) ]
      #points = []
      #points.append((float(rx[0]),float(ry[0]),float(rz[0])))


      for i in range(N):
          p1 =  np.arange(4*i,4*i+4)

          data[p1[0]][0] = rx1[i]
          data[p1[0]][1] = ry1[i]
          data[p1[0]][2] = rz1[i]

          data[p1[1]][0] = rx2[i]
          data[p1[1]][1] = ry2[i]
          data[p1[1]][2] = rz2[i]

          data[p1[2]][0] = rx2[i+1]
          data[p1[2]][1] = ry2[i+1]
          data[p1[2]][2] = rz2[i+1]

          data[p1[3]][0] = rx1[i+1]
          data[p1[3]][1] = ry1[i+1]
          data[p1[3]][2] = rz1[i+1]

          #points.append((float(rx[i+1]),float(ry[i+1]),float(rz[i+1]))) # midline
          #print(points)

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
       Functions ended here. Below is the main program
=========================================================================
'''

'''
=========================================================================
 Creating mobius surface
=========================================================================
'''

#=================
global wd,nfold
wd        = 0.007
thickness = .15*wd
cols = (0,0.976, 0.968,1)  # rgb and facealpha for the mobius surface

'''
=========================================================================
 Name and path of the data file
=========================================================================
'''

branch = 2

if branch==1:

    nfold = 3
    str1 =  '3fold_N'
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/'

    strtau = path+'tau_branch_1.txt'
    tau1 = np.loadtxt(strtau)
    N = 72
    N1=N-1
    h = 1


elif branch==2:

    nfold = 5
    str1 = '5pi_knot_N'
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/'

    strtau = path + 'tau_branch_2.txt'
    tau1 = np.loadtxt(strtau)
    N = 120
    N1= N-1
    h = 1/N
else:

    nfold = 7
    str1 = '7pi_knot_N'
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/'

    strtau = path + 'tau_branch_3.txt'
    tau1 = np.loadtxt(strtau)

#


ln1 = np.size(tau1)

frame_num = 0
for i in range(ln1):

    frame_num= frame_num+1
    tau  =  int(np.around(tau1[i+10]*(10**10),decimals=0))
    #tau  =   np.around(tau1[i+10]*(10**10),decimals=1)

    strb =  path + 'b_' + str1 + str(N)+ '_tau_'+ str(tau)+ '.txt'

    if i==0:
        strb = path + 'b_5pi_knot_N120_tau_119842333032.1598.txt'

    verts = ReadDataFile(strb)
    #   print(type(points))
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


    mySurfaceObject = bpy.data.objects.new("MySurface"+str(i), mesh)

    # add the mesh as an object into the scene with this utility module
    bpy.context.view_layer.active_layer_collection.collection.objects.link(mySurfaceObject)
    mat = bpy.data.materials.new("matBase")
    mat.diffuse_color = cols
    mySurfaceObject.active_material = mat
    mySurfaceObject.material_slots[0].link = 'OBJECT'                               #link material to object not to mesh(data)
    mySurfaceObject.material_slots[0].material = mat



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
    mySurfaceObject.modifiers["Subdivision"].render_levels = 6 #change factor in render
    mySurfaceObject.modifiers["Subdivision"].levels = 6             #change factor in viewport only


    #=========== Shade smooth ======
    bpy.ops.object.editmode_toggle()
    bpy.ops.object.shade_smooth()
    #==============================================================================================

    '''
    ========================= Surface created =====================================
    '''

    '''
    ========================= Creating midline ===================================
    '''
    points = []
    #points = ReadDataFile(userinput['filedata'])
    # Curve C10     - Trivial shape

    colr = (1,0,0,1)#(0.301,0.811,0.498,1)
    #thickness = 0.01
    points = ReadDataFile_r(strb)
    curve1 = CreateCurve(points, thickness, colr)            #create curve



    ob = bpy.data.objects["MySurface"+str(i)]


    #bpy.context.scene.frame_set(0)
    ob.hide_render = True
    ob.keyframe_insert(data_path="hide_render", frame=0)
    #bpy.context.scene.frame_set(i+1)
    ob.hide_render = False
    ob.keyframe_insert(data_path="hide_render", frame=i+1)
    #bpy.context.scene.frame_set(i+2)
    ob.hide_render = True
    ob.keyframe_insert(data_path="hide_render", frame=i+2)



    curve1.hide_render = True
    curve1.keyframe_insert(data_path="hide_render", frame=0)
    #bpy.context.scene.frame_set(i+1)
    curve1.hide_render = False
    curve1.keyframe_insert(data_path="hide_render", frame=i+1)
    #bpy.context.scene.frame_set(i+2)
    curve1.hide_render = True
    curve1.keyframe_insert(data_path="hide_render", frame=i+2)


bpy.context.scene.frame_end = i+1

#objs = bpy.data.objects
#objs.remove(objs["Cube"], do_unlink=True)

#==== Camera location ===
def update_camera(camera, focus_point=mathutils.Vector((0.0, 0.0, 1.0)), distance=0.5):
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
