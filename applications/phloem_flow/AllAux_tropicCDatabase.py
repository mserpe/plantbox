""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np


def toTry():
    
    Qs = np.array([ 120.,  230.,  340.,  450.,  560.,])*1e-6
    Ratiothrshold = np.linspace(0.1,1,6) #np.array([1/10*i for i in range(10)])
    nodeD = np.array([3,4,5,6,7,8])#0
    maxrun = len(Qs) * len(nodeD) * len(Ratiothrshold)


    Qsv, nodeDv,Ratiothrsholdv = np.meshgrid(Qs, nodeD,Ratiothrshold)
    Qsv=Qsv.flatten() 
    nodeDv=nodeDv.flatten()
    Ratiothrsholdv = Ratiothrsholdv.flatten()
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv,  
              'Ratiothrshold':Ratiothrsholdv}
    return dictPara

