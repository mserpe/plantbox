"""plant example"""
import sys
sys.path.append("../../..")
sys.path.append("../../../build")
sys.path.append("../../../src/python_modules")
import numpy as np;
import subprocess;

new_build = False

if False :
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

plant = pb.MappedPlant()

path = "./results/"

plant.readParameters(path + "P0_plant.xml")

time = 28
leaf_res = 30

for p in plant.getOrganRandomParameter(pb.leaf):
  p.lb =  10 # length of leaf stem
  p.lbs =  1 # length of leaf stem
  p.la,  p.lmax = 38.41053981, 38.41053981
  p.areaMax = 54.45388021  # cm2, area reached when length = lmax
  phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
  l = np.array([38.41053981,1 ,1, 0.3, 1, 38.41053981]) #distance from leaf center
  p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
  p.tropismN = 5
  p.tropismS = 0.5
  print(p.dx)
  p.dx = 1
  #p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
  p.createLeafRadialGeometry(phi, l, leaf_res)

for p in plant.getOrganRandomParameter(pb.stem):
  r= 0.07
  p.r = r
  p.lmax = (time-7)*r


# Initialize
plant.initialize()
plant.SetGeometryResolution(8)
plant.SetLeafResolution(leaf_res)

# Simulate
plant.simulate(time, True)

print("This plant has ", plant.getNumberOfNodes(), " nodes")
organs = plant.getOrgans()
print("This plant has ", len(organs), " organs")
if next((o for o in organs if o.getNumberOfNodes() <= 1), False) :
  print("This plant an organ with only one node")
else :
  print("This plant has no organ without nodes")

vp.plot_plant(plant, "organType")