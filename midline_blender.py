#ref for line arguments
#https://docs.blender.org/manual/en/dev/advanced/command_line/arguments.html

'''
#run & render in background:
blender -b sphere.blend -P datatoblender.py filedata=my_data.out basis=0 render=True save-scene=False

#run interactively:
blender sphere.blend -P midline_blender.py filedata=my_data.out basis=0 render=False save-scene=False

'''


import os
import sys
import bpy
import csv
import math
import random
import time
import traceback
import numpy as np


'''
=========================================================================
Returns {} dictionary, with self reserved words.:
GENERIC INPUT
=========================================================================
'''
def ProcessUserInput():
  dict = {}
  na = 0
  for a in sys.argv:
    print("arg" + str(na) + ": " + a)
    spl = a.split("=")
    if len(spl) == 2:
      key = spl[0]
      val = spl[1]
      if val == "None" or val == "":
        val = False
      elif val == "True":
        val = True
      elif val == "False":
        val = False
      dict[key] = val
    na += 1
  return dict




'''
=========================================================================
ReadDatafile
=========================================================================
'''
def ReadDataFile(FilePath):
    points = []
    #with open(FilePath,"r",encoding='utf-8', errors='ignore') as f:
    with open(FilePath,'r') as f:

     # c = csv.reader(f, delimiter=',', skipinitialspace=True)
      for sRow in f:
          a,b,c = sRow.split()                                                           # thats an individual object
          points.append((float(a),float(b),float(c)))

      N = int(len(points)-1)
      xyz = np.array(points)
      print(N)

      bx = np.transpose(np.zeros((N+1)))
      by = np.transpose(np.zeros((N+1)))
      bz = np.transpose(np.zeros((N+1)))

      bx,by,bz = xyz[0:N+1,0],xyz[0:N+1,1],xyz[0:N+1,2]

      '''
      fig1 = plt.figure()
      ax = fig1.add_subplot(111,projection='3d')
      ax.plot(bx,by,bz)
      plt.show()
      '''

      #========================================================
      #               Constructing tangent
      #========================================================
      tx = np.transpose(np.zeros((N+1)))
      ty = np.transpose(np.zeros((N+1)))
      tz = np.transpose(np.zeros((N+1)))

      #== tangent ==
      tx[0:N] = np.multiply(by[0:N],bz[1:N+1]) - np.multiply(bz[0:N],by[1:N+1])
      ty[0:N] = np.multiply(bz[0:N],bx[1:N+1]) - np.multiply(bx[0:N],bz[1:N+1])
      tz[0:N] = np.multiply(bx[0:N],by[1:N+1]) - np.multiply(by[0:N],bx[1:N+1])

      tx[N] = tx[0]#np.insert(tx,N,tx[0,0])
      ty[N] = ty[0]#np.insert(ty,N,ty[0,0])
      tz[N] = tz[0]#np.insert(tz,N,tz[0,0])

      #========================================================

      #========================================================
      #           Constructing midline
      #========================================================
      #== Midline
      rx = np.transpose(np.zeros((N+1)))
      ry = np.transpose(np.zeros((N+1)))
      rz = np.transpose(np.zeros((N+1)))

      for i in range(N):
          rx[i+1] = rx[i] + tx[i]
          ry[i+1] = ry[i] + ty[i]
          rz[i+1] = rz[i] + tz[i]

      rx = rx - np.mean(rx)
      ry = ry - np.mean(ry)
      rz = rz - np.mean(rz)

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
def CreateCurve(points,thickness,color):
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
    mat.diffuse_color = (1,1,0,1)
    curveOBJ.active_material = mat
    curveOBJ.material_slots[0].link = 'OBJECT'                               #link material to object not to mesh(data)
    curveOBJ.material_slots[0].material = mat

    return curveOBJ

'''
==================================================================================================================================================
Program Run
==================================================================================================================================================
'''
try:
    bpy.context.scene.render.use_freestyle = False    #turn off the outline if it is on

    #DrawSphere()   # draw sphere if loading new scene

    userinput = ProcessUserInput()
    points = []
    #points = ReadDataFile(userinput['filedata'])
    # Curve C10     - Trivial shape

    color = (1,0,0)
    thickness = 0.01
    points = ReadDataFile('b_5pi_knot_N120_tau_178.txt')
    curve1 = CreateCurve(points, thickness, color)            #create curve
     
    if 'basis' in userinput and userinput['basis']:
        CreateBasis(points, userinput['basis'],thickness, curve1)
        CreateBasis(points, userinput['basis'],thickness, curve2)

    if 'save-scene' in userinput and userinput['save-scene']:
        bpy.ops.wm.save_as_mainfile(filepath="scene_" + time.strftime("%Y-%m-%d %H:%M") + ".blend")	#save the scene

    if 'render' in userinput and userinput['render']: 							# if we do rendering from python then we quit.
        scene = bpy.context.scene
        scene.render.image_settings.file_format = 'PNG';
        scene.render.antialiasing_samples = '16'
        scene.render.filepath = "out_" + time.strftime("%Y-%m-%d %H:%M") + ".png"
        bpy.ops.render.render( write_still=True )
        bpy.ops.wm.quit_blender()                     							#quit when done



except:
    print(traceback.format_exc())
    raise
#sys.exit()
