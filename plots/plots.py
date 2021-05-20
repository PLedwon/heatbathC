#!/usr/bin/python3
import numpy as np
import glob
import matplotlib.pyplot as plt
import math
import scipy
from scipy.optimize import curve_fit
import pandas as pd

Omega=1.0

if not glob.glob('./*.npz'):

    resultList =glob.glob('../data/*.csv')
    df = pd.read_csv(resultList[0])
    q = df.to_numpy()
    varQ = np.zeros(np.size(q))
    squaredQ = np.zeros(np.size(q))
    aveQ = np.zeros(np.size(q))

    for file in resultList:
        df = pd.read_csv(file)
        q = df.to_numpy()
        squaQ = np.power(q,2.0)
        print(squaQ)
        squaredQ = squaredQ + squaQ
        aveQ = aveQ + q

norm = np.power(len(resultList),-1.0)
varQ = norm*squaredQ - np.power(norm,2.0)*aveQ
plt.plot(varQ)
plt.show()



#    t=data['t']
#
#
##    stdMat = np.zeros((len(t),len(resultList)))
#
#    i=0
#    for file in resultList:
#        results = np.load(file)
##        stdMat[:,i] = np.power(results['Q'],2) - np.average(results[Q])
#        varQ += np.power(results['Q'],2) - np.average(results['Q'])
#        varP += np.power(results['P'],2) - np.average(results['P'])
#        results.close()
#        i+=1
#        print(i)

#    std=np.std(stdMat, axis=1)
#    norm=1.0/(float(len(resultList)))
#    varQ *= norm
#    varP *= norm
#
#########
#    std=0
#    std *= norm
#
#
#    np.savez("./data/data", varQ=varQ,  t=t, std=std, varP=varP )
#
#
#else:
#
#    datafile=glob.glob('/users/stud/ledwon/Documents/data/*.npz')
#    data=np.load(datafile[0])
#    varQ=data['varQ']
#    varP=data['varP']
#    std=data['std']
#    t=data['t']
#
#
#nSaved=5000
#ds=(t[-1]-t[0])/nSaved
#plotTimes=ds*range(1,len(varQ)+1)
#
#print(len(plotTimes),len(varQ))
#

#varQPlot.savefig("./img/varQ.pdf")
#
#
#
#