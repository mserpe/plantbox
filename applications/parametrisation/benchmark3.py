"""Analogy to Shehan's soil core setup"""
import sys;
from pathlib import Path
import os

home = str(Path.home())

CPlantBox_dir = "/4Dirk/CPlantBox"
    
sys.path.append(home +CPlantBox_dir)
sys.path.append(home +CPlantBox_dir+"/src/python_modules")
sys.path.append(home +CPlantBox_dir+"/src")
import numpy as np
#import matplotlib.pyplot as plt
#import vtk_plot as vp
import plantbox as pb
#import pickle as pkl
#import math

import psutil


months=8
times=np.linspace(0,30*months,months+1)
runs = 1 #100 if want to set synthetic objective dataset
rlds_lst = []
dt_ = np.diff(np.asarray(times))

# 72 cm*45 cm size plot
M = 16  # number of plants in rows
N = 7 # number of rows
distp = 3  # distance between the root systems along row[cm]
distr =12  # distance between the rows[cm]
interrow=(M-1)*distp # intra-row spacing
row=(N-1)*distr # row spacing

r, depth, layers = 4.2/2, 160., 32 # Soil core analysis
layerVolume = depth / layers * r * r * np.pi

z_ = np.linspace(0, -depth, layers)  # slices of soil core

soilcolumn = pb.SDF_PlantContainer(r, r, depth, False)#square = false


soilcor_x = interrow/4
soilcor_y = distr
x_= [11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75]                    
y_= [66.0, 66.0, 66.0, 54.0, 54.0, 54.0, 42.0, 42.0, 42.0, 30.0, 30.0, 30.0, 18.0, 18.0, 18.0] 
soilcolumns_ =[pb.Vector3d(x_ij, y_ij,0) for x_ij, y_ij in zip(x_, y_)]
soilcolumns = [pb.SDF_RotateTranslate(soilcolumn, vi) for vi in soilcolumns_]
#for j in reversed(range(5)): #to get right order when plotting
 #   for i in range(3):
  #      x_ = (i+1)*soilcor_x
   #     y_ = j*soilcor_y + distr*1.5
    #    v = pb.Vector3d(x_, y_,0)
    #    colT = pb.SDF_RotateTranslate(soilcolumn, v)
     #   name__ = "col_"+str(i)+"_"+str(j)+".py"
      #  soilcolumns.append(colT)

#gemtry = pb.SDF_Union(soilcolumns)


# creates soil space to stop roots from growing out of the soil
soilSpace = pb.SDF_PlantContainer(500,500, 500, True)

def run_benchmark(name_, filedata): #make a function to be called by DREAM
    if filedata == "":
        name_2 = "/modelparameter - true - data2model/"
        path = name_+ name_2#"/home/rbtlm640/CPlantBox/modelparameter/rootsystem/"
        name = "wheat"#"wheat2"
        inputData = path + name + ".xml"
        fromFile_ = True
    else :
        #name_2 = "/modelparameter/"
        inputData = filedata
        fromFile_ = False
        
    # Initializes N*M root systems
    allRS = []
    #allAna = pb.SegmentAnalyser()#(rs)
    for i in range(0, M):
        for j in range(0, N):
            print("seed could be ", np.random.rand()*100000)
            rs = pb.RootSystem()#np.random.rand()*100000) #random seed set in master file (either run_mcmc or other)
            rs.readParameters(inputData, fromFile = fromFile_)#path + name + ".xml")
            rs.getRootSystemParameter().seedPos = pb.Vector3d(distp * i, distr * j, -3.)  # cm
            #rs.setSeed(0) # to eliminates model stochasticity
            rs.setGeometry(soilSpace) # creates soil space to stop roots from growing out of the soil
            rs.initialize(False)
            for dt in dt_:
                rs.simulate(dt, False)
            #rs.simulate(30*months, False) #why not do it this way?
            allRS.append(rs) #need to save the organism or the weak_ptr in SegmentAnalyser will expire
            if i+ j == 0:
                allAna = pb.SegmentAnalyser(rs) #need to do that first otherwise get errors
            else:
                allAna.addSegments(rs)  # collect all in a segAna object, (slow, adds > 2s)
            #allAna.addSegments(rs) 
            #mid_time = time.time()
            #print("--- %s seconds for plant development---" % (time.time() - start_time))
    # allAna.mapPeriodic(interrow, row) # mapping to periodic domain
    #allAna.write("all_plants.vtp")
    rld=np.zeros([len(soilcolumns)*len(times[1:]), layers])
    maxbound = allAna.getMaxBounds() # check if we have root nodes aboveground
    if maxbound.z > 1.0: #even with boundary, node can be a bit > 0 (current observed max: 0.12 cm) 
        print("highest aboveground node > 1cm :", maxbound)
    
    for k in range(len(soilcolumns)):
        ana = pb.SegmentAnalyser(allAna)  
        #ana.addSegments(allAna)
        ana.crop(soilcolumns[k]) #select soil column
        #ana.write("soilCore_%s.vtp" %(k))
        #ana.pack()
        # print("within benchmark ")
        # process = psutil.Process(os.getpid())
        # print('     curren process',os.getpid(), "resident set size",process.get_memory_info().rss(),
            # "virtual memory",process.get_memory_info().vms(),
            # "memory used by shared libraries",process.get_memory_info().lib())
        for j in range(len(times[1:])):
            ana.filter("creationTime", 0, np.flip(np.asarray(times))[j]) 
            ana.pack()#better efect here seems like?
            distrib = ana.distribution("length", 0., -depth, layers, True)
            #fill RLD matrix
            #print("set soil n*",k, "time", np.flip(np.asarray(times))[j], "at index",(len(times[1:])-1-j)*len(soilcolumns)+k)
            rld[(len(times[1:])-1-j)*len(soilcolumns)+k]= np.array(distrib)/layerVolume
    #print("--- %s seconds for plant development---" % (time.time() - start_time))        
    #raise Exception
    rld = np.transpose(rld)  
    
    # Export container geometry as Paraview Python script
    # gemtry = pb.SDF_Union(soilcolumns)
    # rs.setGeometry(gemtry)
    # rs.write("gemtry.py")
    
    return rld

def forward_model_crr(X):#delayed(forward_process)(X_row,n,extra_par)
    X[0]=np.round(X[0])
    
    template="""<plant>
    <seed name="wheat" subType="0">
        <parameter name="maxTi" value="0"/>
        <!--Maximal number of tillers [1]-->
        <parameter name="delayB" value="0"/>
        <!--Time delay between the basal roots [day]-->
        <parameter name="delayRC" value="1000"/>
        <!--Delay between the root crowns [day]-->
        <parameter name="delaySB" value="1000"/>
        <!--Time delay between the shoot borne roots [day]-->
        <parameter name="firstB" value="0"/>
        <!--Emergence of first basal root [day]-->
        <parameter name="firstSB" value="1000"/>
        <!--First emergence of a shoot borne root [day]-->
        <parameter name="maxB" value="{ppar1:d}"/>
        <!--Maximal number of basal roots [1]-->
        <parameter name="nC" value="0"/>
        <!--Maximal number of roots per root crown [1]-->
        <parameter name="nz" value="0"/>
        <!--Distance between the root crowns along the shoot [cm]-->
        <parameter name="seedPos.x" value="0"/>
        <!--X-Coordinate of seed position [cm]-->
        <parameter name="seedPos.y" value="0"/>
        <!--Y-Coordinate of seed position [cm]-->
        <parameter name="seedPos.z" value="-3"/>
        <!--Z-Coordinate of seed position [cm]-->
        <parameter name="simulationTime" value="240"/>
        <!--Recommended final simulation time  [day]-->
    </seed>
    <root name="taproot" subType="1">
        <parameter name="gf" value="1"/>
        <!--Growth function number [1]-->
        <parameter name="tropismT" value="1"/>
        <!--Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)-->
        <parameter name="a" value="0.050000001"/>
        <!--Root radius [cm]-->
        <parameter name="colorB" value="0.196078"/>
        <!--Root color, blue component [0.-1.]-->
        <parameter name="colorG" value="0.39215699"/>
        <!--Root color, green component [0.-1.]-->
        <parameter name="colorR" value="0.431373"/>
        <!--Root color, red component [0.-1.]-->
        <parameter name="dx" value="0.25"/>
        <!--Axial resolution [cm] (maximal segment size)-->
        <parameter name="la" value="4.2" dev="6.4"/>
        <!--Apical zone [cm]-->
        <parameter name="lb" value="0.8" dev="1.2"/>
        <!--Basal zone [cm]-->
        <parameter name="lmax" value="{rpar3:3.2f}" dev="{rpar4:3.2f}"/>
        <!--Maximal root length [cm]-->
        <parameter name="ln" value="{rpar1:3.2f}" dev="{rpar2:3.2f}"/>
        <!--Inter-lateral distance [cm]-->
        <parameter name="r" value="{rpar5:3.2f}" dev="{rpar6:3.2f}"/>
        <!--Initial growth rate [cm day-1]-->
        <parameter name="rlt" value="1e+09"/>
        <!--Root life time [day]-->
        <parameter name="theta" value="{rpar9:3.2f}" dev="{rpar10:3.2f}"/> 
        <!--Angle between root and parent root [rad]-->
        <parameter name="tropismN" value="{rpar7:3.2f}"/> 
        <!--Number of trials of root tropism-->
        <parameter name="tropismS" value="{rpar8:3.2f}"/>
        <!--Mean value of expected change of root tropism [1/cm]-->
        <parameter name="successor" type="2" percentage="1"/>
        <!--Sub type of lateral roots-->
    </root>
    <root name="1storderlateral" subType="2">
        <parameter name="gf" value="1"/>
        <!--Growth function number [1]-->
        <parameter name="tropismT" value="1"/>
        <!--Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)-->
        <parameter name="a" value="0.029999999"/>
        <!--Root radius [cm]-->
        <parameter name="colorB" value="0.156863"/>
        <!--Root color, blue component [0.-1.]-->
        <parameter name="colorG" value="0.49019599"/>
        <!--Root color, green component [0.-1.]-->
        <parameter name="colorR" value="0.54901999"/>
        <!--Root color, red component [0.-1.]-->
        <parameter name="dx" value="0.25"/>
        <!--Axial resolution [cm] (maximal segment size)-->
        <parameter name="la" value="{rpar11:3.2f}" dev="{rpar12:3.2f}"/>
        <!--Apical zone [cm]-->
        <parameter name="lb" value="0.80" dev="1."/>
        <!--Basal zone [cm]-->
        <parameter name="lmax" value="{rpar15:3.2f}" dev="{rpar16:3.2f}"/>
        <!--Maximal root length [cm]-->
        <parameter name="ln" value="{rpar13:3.2f}" dev="{rpar14:3.2f}"/>
        <!--Inter-lateral distance [cm]-->
        <parameter name="r" value="0.40" dev="0.12"/>
        <!--Initial growth rate [cm day-1]-->
        <parameter name="rlt" value="1e+09"/>
        <!--Root life time [day]-->
        <parameter name="theta" value="1.2" dev="0.4000"/>
        <!--Angle between root and parent root [rad]-->
        <parameter name="tropismN" value="1"/>
        <!--Number of trials of root tropism-->
        <parameter name="tropismS" value="0.4000"/>
        <!--Mean value of expected change of root tropism [1/cm]-->
        <parameter name="successor" type="3" percentage="1"/>
        <!--Sub type of lateral roots-->
    </root>
    <root name="2ndorderlateral" subType="3">
        <parameter name="gf" value="1"/>
        <!--Growth function number [1]-->
        <parameter name="tropismT" value="1"/>
        <!--Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)-->
        <parameter name="a" value="0.02"/>
        <!--Root radius [cm]-->
        <parameter name="colorB" value="0.117647"/>
        <!--Root color, blue component [0.-1.]-->
        <parameter name="colorG" value="0.58823502"/>
        <!--Root color, green component [0.-1.]-->
        <parameter name="colorR" value="0.627451"/>
        <!--Root color, red component [0.-1.]-->
        <parameter name="dx" value="0.25"/>
        <!--Axial resolution [cm] (maximal segment size)-->
        <parameter name="la" value="2.2" dev="0.4"/>
        <!--Apical zone [cm]-->
        <parameter name="lb" value="0"/>
        <!--Basal zone [cm]-->
        <parameter name="lmax" value="2.18" dev="0.56"/>
        <!--Maximal root length [cm]-->
        <parameter name="ln" value="0"/>
        <!--Inter-lateral distance [cm]-->
        <parameter name="r" value="1.0" dev="0.2"/>
        <!--Initial growth rate [cm day-1]-->
        <parameter name="rlt" value="1e+09"/>
        <!--Root life time [day]-->
        <parameter name="theta" value="1.12" dev="0.40"/>
        <!--Angle between root and parent root [rad]-->
        <parameter name="tropismN" value="0.1"/>
        <!--Number of trials of root tropism-->
        <parameter name="tropismS" value="0.600"/>
        <!--Mean value of expected change of root tropism [1/cm]-->
    </root>
    </plant>"""
    context= {"ppar1":X[0].astype('int'),#NB
               "rpar1":X[1],#ln0 ln0s
                  "rpar2":X[2],
                  "rpar3":X[3],#lmax0, lmax0s
                  "rpar4":X[4],
                  "rpar5":X[5],#r0, r0s
                  "rpar6":X[6],
                  "rpar7":X[7],#tr0, tr0s 
                  "rpar8":X[8],
                  "rpar9":X[9],#theta0, theta0s
                  "rpar10":X[10],
                  "rpar11":X[11],#la1, la1s
                  "rpar12":X[12],
                  "rpar13":X[13],#ln1, ln1s
                  "rpar14":X[14],
                  "rpar15":X[15],#maxl1, maxl1s
                  "rpar16":X[16]
                  }
    name = ""
    fx_ = run_benchmark(name, template.format(**context))
    fx= fx_.flatten()
    return fx

if __name__ == '__main__': #for testing
    import time
    start_time = time.time()
    rlds = run_benchmark(".", "")
    np.savetxt("."+'/wheat_core_matrix.txt', rlds)#for testing   
    print("--- %s seconds, end benchmark ---" % (time.time() - start_time))
    np.set_printoptions(threshold=sys.maxsize)
    
    
    import matplotlib.pyplot as plt
    import vtk_plot as vp
    import pickle as pkl
    import math
    # plotting RLDs
    rlds_mean = rlds
    
    #re-organise rld for plotting:
    rlds_mean = np.transpose(rlds_mean)
        
        
    pt_idx=[(0,0),(0,1),(0,2),
    (1,0),(1,1),(1,2),
    (2,0),(2,1),(2,2),
    (3,0),(3,1),(3,2),
    (4,0),(4,1),(4,2),]
    #legend_lst=[]
    #for i in range(len(times)-1):
    #    legend_lst.append(str(int(times[i+1])))
    legend_lst = [str(int(t_i)) for t_i in times[1:]]    
    fig, axes = plt.subplots(nrows = 5, ncols = int(len(soilcolumns)/5), sharex=True, sharey=True,figsize = (8,16))
    
    
    for k in range(len(soilcolumns)):
        axes[pt_idx[k]].set_title('Soil core'+' ' +str(k+1))
        for j in range(len(times[1:])):
            #print("for soil n*",k,"time",legend_lst[j],"index", len(soilcolumns)*j+k)
            #print(np.array(rlds_mean[len(soilcolumns)*j+k]))
            axes[pt_idx[k]].plot(np.array(rlds_mean[len(soilcolumns)*j+k]), z_)
            axes[pt_idx[k]].set_xlim(0,5)
            
    plt.setp(axes[-1, :], xlabel='RLD $(cm/cm^3)$')
    plt.setp(axes[:, 0], ylabel='Depth $(cm)$')
    plt.legend(np.asarray(legend_lst),loc="lower center", bbox_to_anchor=(-0.8, -0.5), ncol= 8)
    fig.subplots_adjust()
    #plt.savefig("rld_plot.png",dpi=300, bbox_inches='tight')
    plt.show()
    