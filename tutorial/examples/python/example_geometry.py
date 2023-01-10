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

path = "../../../modelparameter/plant/"

plant.readParameters(path + "wheat_withStemLeaf.xml")

# Initialize
plant.initialize()

# Simulate
plant.simulate(1, True)

print("test")

vtkplant = vtk.vtkPolyData()
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
geom = np.reshape(geom, (int(geom.shape[0]/3), 3))
print(geom)
points.SetData(numpy_to_vtk(geom))
cell_data = np.array(plant.GetGeometryIndices())
cell_data = np.reshape(cell_data, (cell_data.shape[0]//3, 3))
cells = numpy_to_vtk(cell_data)
pd.SetPoints(points)
pd.SetPolys(cells)

print("Extracted the Plant")

print("test")

print(vtkplant.GetPoints())
