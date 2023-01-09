"""plant example"""
import sys
sys.path.append("../../..")
sys.path.append("../../../build")
sys.path.append("../../../src/python_modules")
import numpy as np;

import plantbox as pb
import vtk_plot as vp

plant = pb.MappedPlant()

path = "../../../modelparameter/plant/"

plant.readParameters("wheat_withStemLeaf.xml")

# Initialize
plant.initialize()

# Simulate
plant.simulate(30, True)

vtkplant = vp.convert_vis_to_vtk(plant)


