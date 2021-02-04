import sys;
sys.path.append("../../..")
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from pyevtk.hl import gridToVTK
import plantbox as rb

rg_ = [0.1, 0.5, 1, 1.5]
for jj in range(0,len(rg_)):

    #
    # Root system
    #
    rs = rb.RootSystem()

    path = "../../../modelparameter/rootsystem/"
    name = "straight_root"  
    rs.readParameters(path + name + ".xml")

    rg = rg_[jj]
    for p in rs.getRootRandomParameter():
        p.r = rg  # growth rate

    #set geometry 
    width = 5  # cm
    depth = 25    
    soilcore = rb.SDF_PlantContainer(width, width, depth, True)
    rs.setGeometry(soilcore)  
    rs.setSeed(0)

    rs.initialize()

    for ii in range(0,10): 
        simtime = 1
        rs.simulate(simtime, True);
    rs.write("vtp/faba_citrate/rg"+str(rg) +"_day"+("{:02d}".format(ii+1))+".vtp")
    #
    # Grid parameter
    #
    nodes = np.array([np.array(n) for n in rs.getNodes()])
    np.save("nodes/faba_citrate/rg"+str(rg) + "_day"+str(ii+1), nodes)
    xres = 0.05
    yres = 0.05
    zres = 0.05
    nx = int(width / xres);
    ny = int(width / yres);
    nz = int(depth / zres);
    print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)

    #
    # Model parameter
    #
    model = rb.ExudationModel2(width, width, depth, nx, ny, nz, rs)
    model.Q = 18.4  # µg/d/tip
    model.Dl = 0.171  # cm2/d
    model.theta = 0.3	 #-
    model.R = 16.7  # -
    model.k = 1.42  # d-1
    model.l = 5  # cm (for line source only)

    #
    # Numerical parameter
    #
    model.type = rb.IntegrationType.mls;  # mps, mps_straight, mls
    model.n0 = 20  # integration points per cm
    model.thresh13 = 1.e-15;  # threshold to neglect diffusing g (eqn 13)
    model.calc13 = True;  # turns Eqn 13  on (True) and off (False)
    model.observationRadius = 5;  # limits computational domain around roots [cm]

    t = time.time()
    model.makeVoxelLists() 
    C = model.calculate(ii+1)
    elapsed = time.time() - t
    print("Computation took", elapsed, "s")

    #
    # post processing...
    #
    C_ = np.zeros((nx,ny,nz))
    C = np.reshape(C, (nz, ny, nx)) 
    for i in range(0,np.shape(C)[0]):
        for j in range(0,np.shape(C)[1]):
            for k in range(0,np.shape(C)[2]):
                C_[k,j,i] = C[i,j,k]
    del C
    C = C_

    np.save("concentration/faba_citrate/rg"+str(rg) +"_day"+("{:02d}".format(ii+1)), C)

    X = np.linspace(-width / 2, width / 2, nx)
    Y = np.linspace(-width / 2, width / 2, ny)
    Z = np.linspace(-depth, 0, nz)

    gridToVTK("exud/faba_citrate/./Exudates_rg"+str(rg) + "_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":C})