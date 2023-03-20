"""plant example"""
import sys
sys.path.append("../../..")
sys.path.append("../../../build")
sys.path.append("../../../src/python_modules")
import numpy as np;
import subprocess;

new_build = False

if new_build :
  import os
  cwd = os.getcwd()
  print("We are in ", cwd)
  # switch to the directory where the plantbox library is located
  os.chdir("../../../build")
  subprocess.call("make -j8 plantbox", shell=True)
  os.chdir(cwd)

import plantbox as pb
import vtk_plot as vp
from vtk.util.numpy_support import *
import vtk


path = "./results/"

# set parameter file range
P_files = [2,3]
P_names = ["P"+str(i)+"_plant" for i in P_files]

time = 28
leaf_res = 30

for versuch in range (0,100) :
  for k,name in enumerate(P_names) :
    plant = pb.MappedPlant()
    plant.readParameters(path + name + ".xml")
    if name == "P0_plant":
      for p in plant.getOrganRandomParameter(pb.leaf):
        p.la,  p.lmax = 38.41053981, 38.41053981
        #p.theta = 0.2 # 0.2
        p.theta = 0.01
        p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
        p.areaMax = 54.45388021  # cm2, area reached when length = lmax
        phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi    
        l = np.array([38.41053981,1 ,1, 0.3, 1, 38.41053981]) #distance from leaf center
        p.leafGeometryPhi = phi
        p.leafGeometryX = l
        #p.tropismN = 5
        p.tropismS = 0.08
        p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
        print("test")
        p.createLeafRadialGeometry(phi,l,leaf_res)
      for p in plant.getOrganRandomParameter(pb.stem):
        p.r = 0.758517633
        p.lb = 4
        p.delayLat = 4
        p.lmax = (time-7)*p.r
        p.dx = 0.1
    elif name == "P1_plant":
      for p in plant.getOrganRandomParameter(pb.leaf):
        p.lb =  0 # length of leaf stem
        p.la,  p.lmax = 42.60617256, 42.60617256
        p.areaMax = 66.69532685  # cm2, area reached when length = lmax
        phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
        l = np.array([42.60617256,1 ,1, 0.3, 1, 42.60617256]) #distance from leaf center
        p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
        #p.tropismN = 5
        p.theta = 0.1
        p.tropismS = 0.08
        p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
        p.createLeafRadialGeometry(phi, l, leaf_res)
      for p in plant.getOrganRandomParameter(pb.stem):
        p.r= 0.91546738
        p.lb = 6
        p.ln = 3
        p.delayLat = 4
        p.lmax = (time-7)*p.r  
    elif name == "P2_plant":
      for p in plant.getOrganRandomParameter(pb.leaf):
        p.lb =  0 # length of leaf stem
        p.la,  p.lmax = 52.23664394, 52.23664394
        p.areaMax = 80.68274258  # cm2, area reached when length = lmax
        phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
        l = np.array([52.23664394,1 ,1, 0.3, 1, 52.23664394]) #distance from leaf center
        p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
        #p.tropismN = 5
        p.tropismS = 0.04
        p.theta = 0.1
        p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
        p.createLeafRadialGeometry(phi, l, leaf_res)
      for p in plant.getOrganRandomParameter(pb.stem):
        p.r= 1.000613891
        p.lb = 6
        p.ln = 3
        p.delayLat = 4
        p.lmax = (time-7)*p.r 
    elif name == "P3_plant":
      for p in plant.getOrganRandomParameter(pb.leaf):
        p.lb =  0 # length of leaf stem
        p.la,  p.lmax = 49.12433414, 49.12433414
        p.areaMax = 71.95670914  # cm2, area reached when length = lmax
        phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
        l = np.array([49.12433414,1 ,1, 0.3, 1, 49.12433414]) #distance from leaf center
        p.theta = 0.1
        p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
        p.tropismN = 5
        p.tropismS = 0.04
        p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
        p.createLeafRadialGeometry(phi, l, leaf_res)
      for p in plant.getOrganRandomParameter(pb.stem):
        p.r= 1.128705967
        p.lb = 6
        p.ln = 3
        p.delayLat = 4
        p.lmax = (time-7)*p.r  
    else :
      exit(-1)
    #
    #plant.readParameters("results/P0_plant.xml")



    # Initialize
    plant.initialize()
    plant.SetGeometryResolution(8)
    plant.SetLeafResolution(leaf_res)

    # Simulate
    plant.simulate(time, True)
    #vp.plot_plant(plant,)

    print("This plant has ", plant.getNumberOfNodes(), " nodes")
    organs = plant.getOrgans()
    print("This plant has ", len(organs), " organs")
    if next((o for o in organs if o.getNumberOfNodes() <= 1), False) :
      print("This plant an organ with only one node")
    else :
      print("This plant has no organ without nodes")

    for organ_type in [3,4] :
      print("Created the plant")
      print("Created Polydata")
      pd = vtk.vtkPolyData()
      # clear the data
      pd.Reset()
      print("Created Points")
      points = vtk.vtkPoints()
      # clear the points
      points.Reset()
      print("Computing Geometry")
      #plant.ComputeGeometryForOrganType(pb.leaf)
      #plant.ComputeGeometry()
      print("Extracting Data")
      plant.SetComputeMidlineInLeaf(organ_type != 4)
      print("Extracting Data from ", organ_type)
      plant.ComputeGeometryForOrganType(organ_type)
      geom = np.array(plant.GetGeometry())
      print(geom.shape)

      print("Extracted Data")

      print("Creating VTK Data from node ids")
      nodeids = np.array(plant.GetGeometryNodeIds())
      print(nodeids.shape)
      nodeids = numpy_to_vtk(nodeids, deep=True)
      nodeids.SetName("nodeids")
      pd.GetPointData().AddArray(nodeids)
      texcoords = np.array(plant.GetGeometryTextureCoordinates())
      print(texcoords.shape)
      texcoords = np.reshape(texcoords, (texcoords.shape[0]//2, 2))
      texcoords = numpy_to_vtk(texcoords, deep=True)
      texcoords.SetName("tcoords")
      pd.GetPointData().AddArray(texcoords)
      normals = np.array(plant.GetGeometryNormals())
      print(normals.shape)
      normals = np.reshape(normals, (normals.shape[0]//3, 3))
      normals = numpy_to_vtk(normals, deep=True)
      normals.SetName("normals")
      pd.GetPointData().AddArray(normals)
      points = vtk.vtkPoints()
      points.SetNumberOfPoints(geom.shape[0]//3)
      for i in range(geom.shape[0]//3) :
        points.SetPoint(i, geom[i*3], geom[i*3+1], geom[i*3+2])
      cell_data = np.array(plant.GetGeometryIndices())
      cell_data = np.reshape(cell_data, (cell_data.shape[0]//3, 3))
      cells = vtk.vtkCellArray()
      for i in range(cell_data.shape[0]) :
        cells.InsertNextCell(3)
        cells.InsertCellPoint(cell_data[i, 0])
        cells.InsertCellPoint(cell_data[i, 1])
        cells.InsertCellPoint(cell_data[i, 2])
      pd.GetPointData().SetTCoords(texcoords)
      pd.GetPointData().SetScalars(nodeids)
      pd.GetPointData().SetNormals(normals)
      pd.SetPoints(points)
      pd.SetPolys(cells)
      # Calculate surface tangents
      tangents = vtk.vtkPolyDataTangents()
      tangents.SetInputData(pd)
      tangents.Update()

      print(pd.GetPoints())
      writer = vtk.vtkXMLPolyDataWriter()
      writer.SetFileName("plants/"+name + "_" + str(time) + "_" + str(organ_type)+"_"+str(versuch) + ".vtp")
      writer.SetDataModeToAscii()
      writer.SetInputData(tangents.GetOutput())
      writer.Write()
    del plant
