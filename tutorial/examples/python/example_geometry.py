"""plant example"""
import sys
sys.path.append("../../..")
sys.path.append("../../../build")
sys.path.append("../../../src/python_modules")
import numpy as np;

import plantbox as pb
import vtk_plot as vp
from vtk.util.numpy_support import *
import vtk

plant = pb.MappedPlant()

path = "./results/"

plant.readParameters(path + "P0_plant.xml")

# a function that filters an array by


# Initialize
plant.initialize()
plant.SetGeometryResolution(8)
plant.SetLeafResolution(20)

# Simulate
plant.simulate(30, True)

print("This plant has ", plant.getNumberOfNodes(), " nodes")
organs = plant.getOrgans()
print("This plant has ", len(organs), " organs")
if next((o for o in organs if o.getNumberOfNodes() <= 1), False) :
  print("This plant an organ with only one node")
else :
  print("This plant has no organ without nodes")

#plot cpb plant
plant.write("test.vtp")

print("test")

print("Created the plant")

print("Created Polydata")
pd = vtk.vtkPolyData()
print("Created Points")
points = vtk.vtkPoints()
print("Computing Geometry")
plant.ComputeGeometry()
print("Extracting Data")
geom = np.array(plant.GetGeometry())
print(geom.shape)

print("Extracted Data")

print("Creating VTK Data from node ids")
nodeids = np.array(plant.GetGeometryNodeIds())
print(nodeids.shape)
nodeids = numpy_to_vtk(nodeids, deep=True)
nodeids.SetName("nodeids")
print("Adding VTK Data to Polydata")
pd.GetPointData().AddArray(nodeids)

texcoords = np.array(plant.GetGeometryTextureCoordinates())
print(texcoords.shape)
texcoords = numpy_to_vtk(texcoords, deep=True)
texcoords.SetName("texcoords")
pd.GetPointData().AddArray(texcoords)

normals = np.array(plant.GetGeometryNormals())
print(normals.shape)
normals = numpy_to_vtk(normals, deep=True)
normals.SetName("normals")
pd.GetPointData().AddArray(normals)

print("Iterating through points to create points form geometry")
points = vtk.vtkPoints()
print("Created point array")
print(geom.shape[0]//3, " points")
points.SetNumberOfPoints(geom.shape[0]//3)
print("Setting ", geom.shape[0]//3, " points")
for i in range(geom.shape[0]//3) :
  points.SetPoint(i, geom[i*3], geom[i*3+1], geom[i*3+2])

print("Getting the cells from the vis")
cell_data = np.array(plant.GetGeometryIndices())
cell_data = np.reshape(cell_data, (cell_data.shape[0]//3, 3))

print("doing the same with the cells")
cells = vtk.vtkCellArray()
for i in range(cell_data.shape[0]) :
  cells.InsertNextCell(3)
  cells.InsertCellPoint(cell_data[i, 0])
  cells.InsertCellPoint(cell_data[i, 1])
  cells.InsertCellPoint(cell_data[i, 2])


print("Adding points and cells to polydata")
pd.SetPoints(points)
pd.SetPolys(cells)

print("Extracted the Plant")

print("test")

print(pd.GetPoints())
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName("test_geom.vtp")
writer.SetDataModeToAscii()
writer.SetInputData(pd)
writer.Write()

print("Written the Plant")

print("Sanity check")
print(pd.GetNumberOfPoints(), " points and ", pd.GetNumberOfCells(), " cells")
print("Maximum cell id: ", cell_data.max())
print("Minimum cell id: ", cell_data.min())
