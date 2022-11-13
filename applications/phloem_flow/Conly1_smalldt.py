""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

from datetime import datetime, timedelta
from masterI import runSim
from ConlyDatabase import toTry

main_dir=os.environ['PWD']#dir of the file
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
results_dir = main_dir +"/results"+directoryN


if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    print(len(test))
    for item in test:
        if item != '.ipynb_checkpoints':
            os.remove(item)
    test = os.listdir(results_dir)
    print(len(test))
    if len(test) >1:
        print("dir not empty")
        raise Exception
        
isCluster = (os.environ['HOME'] == '/home/m.giraud')
maxcore = 8
if isCluster:
    maxcore = 256
parallelizer = Parallel(n_jobs=maxcore - 1)



params = toTry()
Qsv=params['Qsv']
MulimSucv=params['MulimSucv']
nodeDv=params['nodeDv']
totrun = 0
maxrun = len(Qsv)

n_jobs = min((maxrun - totrun),
                 maxcore - 1)
print("leftToDo:", maxrun - totrun)
tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = 0, thresholdAux = 0, 
                         RatiothresholdAux =10,
       Qmax_ = Qsv[i+totrun], thresholdSuc = MulimSucv[i+totrun], 
       useCWGr = True, UseRatiothresholdAux = False,
                         nodeD = nodeDv[i+totrun], thread = i,
       activeAtThreshold_auxin = False, activeAtThreshold_suc = True,
                         testTime=7, dtBefore = 1/24, dtAfter= 1/(24*60))
                    for i in range(n_jobs))
parallelizer(tasks_iterator)
totrun += n_jobs