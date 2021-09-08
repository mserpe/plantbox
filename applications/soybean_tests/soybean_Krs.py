""" 
soil uptake fraction of a root system (soil is in hydrostatic equilibrium) 
"""
import sys; sys.path.append("../../python/modules/"); sys.path.append("../../../CPlantBox/"); 
from cmath import isnan
sys.path.append("../../../CPlantBox/src/python_modules")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
# artificial shoot
kr0 = np.array([[0, 1.e-12], [1e20, 1.e-12]])		# kr is very small
kz0 = np.array([[0, 0.356832], [1e20, 0.356832]])	# kx is equal to tap root

# tap root
kr1 = np.array([[0, 1.14048e-3], [2, 1.08864e-3], [4, 1.0368e-3], [6, 9.8486e-4], [8, 9.3312e-4], [10, 8.8992e-4], [12, 8.47584e-4], [14, 8.06112e-4], [16, 7.67232e-4], [18, 7.3008e-4], [20, 6.9552e-4], [22, 6.61824e-4], [24, 6.29856e-4], [26, 5.99616e-4], [28, 5.7024e-4], [30, 5.42592e-4], [32, 5.16672e-4], [1e20, 5.16672e-4]])
kz1 = np.array([[0, 0.067392], [2, 0.074736], [4, 0.082944], [6, 0.092448], [8, 0.101952], [10, 0.113184], [12, 0.126144], [14, 0.139968], [16, 0.154656], [18, 0.171936], [20, 0.190944], [22, 0.21168], [24, 0.235008], [26, 0.260928], [28, 0.28944], [30, 0.321408], [32, 0.356832], [1e20, 0.356832]])

# first order lateral
kr2 = np.array([[0, 4.11264e-3], [1, 3.888e-3], [2, 3.672e-3], [3, 3.47328e-3], [4, 3.2832e-3], [5, 3.10176e-3], [6, 2.92896e-3], [7, 2.77344e-3], [8, 2.61792e-3], [9, 2.47968e-3], [10, 2.34144e-3], [11, 2.21184e-3], [12, 2.09952e-3], [13, 1.97856e-3], [14, 1.86624e-3], [15, 1.76256e-3], [16, 1.66752e-3], [17, 1.58112e-3], [1e20, 1.58112e-3]])
kz2 = np.array([[0, 4.06944e-4], [1, 5.00256e-4], [2, 6.15168e-4], [3, 7.56e-4], [4, 9.3312e-4], [5, 1.14048e-3], [6, 1.40832e-3], [7, 1.728e-3], [8, 2.12544e-3], [9, 2.60928e-3], [10, 3.21408e-3], [11, 3.94848e-3], [12, 4.85568e-3], [13, 5.97024e-3], [14, 7.344e-3], [15, 8.9856e-3], [16, 0.0110592], [17, 0.0136512], [1e20, 0.0136512]])

# second order lateral
kr3 = np.array([[0, 4.11264e-3], [1, 3.888e-3], [2, 3.672e-3], [3, 3.47328e-3], [4, 3.2832e-3], [5, 3.10176e-3], [6, 2.92896e-3], [7, 2.77344e-3], [8, 2.61792e-3], [9, 2.47968e-3], [10, 2.34144e-3], [11, 2.21184e-3], [12, 2.09952e-3], [13, 1.97856e-3], [14, 1.86624e-3], [15, 1.76256e-3], [16, 1.66752e-3], [17, 1.58112e-3], [1e20, 1.58112e-3]])
kz3 = np.array([[0, 4.06944e-4], [1, 5.00256e-4], [2, 6.15168e-4], [3, 7.56e-4], [4, 9.3312e-4], [5, 1.14048e-3], [6, 1.40832e-3], [7, 1.728e-3], [8, 2.12544e-3], [9, 2.60928e-3], [10, 3.21408e-3], [11, 3.94848e-3], [12, 4.85568e-3], [13, 5.97024e-3], [14, 7.344e-3], [15, 8.9856e-3], [16, 0.0110592], [17, 0.0136512], [1e20, 0.0136512]])

path = "../../../CPlantBox//modelparameter/rootsystem/"
name = "Moraesetal_2020"  # ""

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")

# Create and set geometry
rs.setMinDx(1.e-3)
x0 = pb.Vector3d(0., 0., -1.)
nx = pb.Vector3d(1., 0., -1.)
ny = pb.Vector3d(0., 1., -1.)
soil_layer = pb.SDF_HalfPlane(x0, nx, ny)  # there was bug, with updated CPlantBox
rs.setGeometry(soil_layer)

rs.setSeed(2)
rs.initialize()

simtime = 154  # days
dt = 1.
N = round(simtime / dt)  # steps

# Plot some length over time
stype = "length"
v_, v1_, v2_, v3_, v4_ = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
segment_ = np.zeros(N)
for i in range(0, N):
    rs.simulate(dt)
    t = np.array(rs.getParameter("type"))
    v = np.array(rs.getParameter(stype))
    segment = np.array(rs.getNumberOfSegments())
    v_[i] = np.sum(v)
    v1_[i] = np.sum(v[t == 1])
    v2_[i] = np.sum(v[t == 2])
    v3_[i] = np.sum(v[t == 3])
    v4_[i] = np.sum(v[t == 4])
    segment_[i] = np.sum(segment)
print("Total number of segments: ", segment_[-1])

""" root system """
rs = pb.MappedRootSystem()
rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, 1.e6))  # not allowed to grow upwards out of soil
p_s = np.linspace(-500, -200, 3001)  #  -200.*np.ones((2001, 1))   # 3 meter down, from -200 to -500, resolution in mm
soil_index = lambda x, y, z: int(-10 * z)  # maps to p_s (hydrostatic equilibirum)
rs.setSoilGrid(soil_index)

path = "../../../CPlantBox//modelparameter/rootsystem/"
name = "Moraesetal_2020"  # "Glycine_max"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
rs.getRootSystemParameter().seedPos.z = -3
rs.setSeed(2)  # random
rs.initialize()

kr4 = kr1  # basal
kz4 = kz1
kr5 = kr1  # shoot borne
kz5 = kz1

krs_, suf_, jc_ = [], [],  []
""" numerical solution of transpiration -1 cm3/day"""
for j in range(1, simtime):

	rs.simulate(1, False)
	
	""" set up xylem parameters """
	r = XylemFluxPython(rs)
	r.setKrTables([kr0[:, 1], kr1[:, 1], kr2[:, 1], kr3[:, 1], kr4[:, 1], kr5[:, 1]], [kr0[:, 0], kr1[:, 0], kr2[:, 0], kr3[:, 0], kr4[:, 0], kr5[:, 0]])
	r.setKxTables([kz0[:, 1], kz1[:, 1], kz2[:, 1], kz3[:, 1], kz4[:, 1], kz5[:, 1]], [kz0[:, 0], kz1[:, 0], kz2[:, 0], kz3[:, 0], kz4[:, 0], kz5[:, 0]])

	suf = r.get_suf(j)
	print("Sum of SUF", np.sum(suf), "from", np.min(suf), "to", np.max(suf), "summed positive", np.sum(suf[suf >= 0]))
	
	krs, jc = r.get_krs(j)
	print("Krs: ", krs)
	print("time: ", j)
	
	krs_.append(krs)
	suf_.append(suf)
	jc_.append(jc)

fig, ax1 = plt.subplots()
t_ = np.linspace(dt, N * dt, N)
ax1.plot(t_, v_, t_, v1_, t_, v2_, t_, v3_, t_, v4_)
ax1.set_xlabel("Time [days]")
ax1.set_ylabel(stype + " [cm]")
#ax1.set_ylim(0, 8000)
ax1.legend(["total", "tap root", "lateral", "2. order lateral", "basal root"], loc = "upper left")

ax1b = ax1.twinx()
time = np.linspace(0, simtime - 1, simtime - 1)
ax1b.plot(time, krs_, "k-")
ax1b.set_ylabel("Global root system conductance $[cm^2 d^{-1}]$")
ax1b.legend(["$K_{rs}$"], loc = "lower right")
plt.savefig("results/" + name + "/" + name + "_Krs.pdf", dpi = 300, bbox_inches='tight')
plt.show()
