import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import plantbox as pb
import vtk_plot as vp # for quick 3d vizualisations
import matplotlib.pyplot as plt # for 2d plots
import numpy as np

# Create instance describing a root system
plant = pb.Plant()

# Open plant and root parameter from a file
path = "../../../modelparameter/plant/"
name = "0"#"UQpresentation"
plant.readParameters(path + name + ".xml")

# define 2D shape
for p in plant.getOrganRandomParameter(pb.leaf):
    p.a = 0.05
    p.a_s = 0
    if (p.subType > 0): 
        #print(p.subType, "radius", p.a, "lmax", p.lmax, p.ln, p.lb, p.successor, p.nob())        
        if (p.subType > 2): 
            #print(p)
            
            p.la, p.lb, p.lmax, p.ln, = 3.5, 1., 7.5, 3  
            p.areaMax = 10  # cm2
            phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
            l= np.array([3, 2.2, 1.7, 2, 3.5])
            N = 101  # N is rather high for testing
            p.createLeafRadialGeometry(phi, l, N)
            p.parametrisationType = 0
            p.tropismT = 1
            p.tropismN = 5
            p.tropismS = 0.1
   
        else:
            p.a = p.a * 3

# Simulate
plant.initialize() 
plant.simulate(30) # [days]

# Export
plant.write("first_example.vtp") # for Paraview
plant.write("first_example.rsml") # e.g. gui/viewer/rsml_viewer.py

# Visualize
_ = vp.plot_plant(plant, "subType") # Plot, using vtk (e.g. "subType")


print("2D leaf shape of a full grown leaf")
lorg = plant.getOrgans(pb.leaf)[1]
lrp = lorg.getLeafRandomParameter()    
leafRadial = (lrp.parametrisationType == 0)
if leafRadial:
    N = len(lrp.leafGeometry)
    yy = np.linspace(0, lorg.leafLength(), N)
    geom_x, geom_y = [],[]
    for i, x in enumerate(lrp.leafGeometry):
        geom_x.extend(x)
        geom_y.extend([yy[i]] * len(x))
    geom_x = np.array(geom_x)
    geom_y = np.array(geom_y)        
    a  = lorg.leafArea() / lorg.leafLength() # scale radius
    plt.plot(geom_x * a, geom_y, "g*")
    plt.plot(-geom_x * a, geom_y, "g*")

else:
    geom_x_a =  np.array([0])
    geom_x_b = np.array([ x[-1] for x in lrp.leafGeometry]) #normalized x value along length
    geom_x = np.concatenate((geom_x_a,geom_x_b))
    geom_y_a = np.array([0])
    geom_y_b =np.linspace(lrp.lb, lorg.leafLength()+lrp.lb, len(geom_x_b))
    geom_y = np.concatenate((geom_y_a,geom_y_b))
    a  = lorg.leafArea() / lorg.leafLength() # scale radius
    plt.plot(geom_x * a, geom_y, "g-*")
    plt.plot(-geom_x * a, geom_y, "g-*")
plt.ylim([0, lrp.lmax+1])
plt.xlim([-a-1, a+1])
plt.axis('scaled')
plt.show()