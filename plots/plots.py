#!/usr/bin/python3
import numpy as np
import glob
import matplotlib.pyplot as plt
import math
import scipy
from scipy.optimize import curve_fit
from random import randrange
import pandas as pd

Omega=1.0

if not glob.glob('./*.npz'):

    resultList =glob.glob('../csvData/*.csv')
    df = pd.read_csv(resultList[0])
    q = df.to_numpy()
    varQ = np.zeros(len(q)-1)
    squaredQ = np.zeros(len(q)-1)
    aveQ = np.zeros(len(q)-1)

    for file in resultList:
        df = pd.read_csv(file)
        q = df.to_numpy()
        print(len(q))
        squaQ = np.power(q,2)
        for i in range(len(squaredQ)):
            squaredQ[i] +=q[i]**2
            aveQ[i] += q[i]

#trajectories = np.zeros(len(q)-1,6)
#    for i in range(0,6):
#        random_index = randrange(len(resultList))
#        df = pd.read_csv(resultList[0])
#        trajectories[:,i]= df.to_numpy()
   

 


norm = np.power(len(resultList),-1.0)
for i in range(len(varQ)):
    varQ[i] = norm*squaredQ[i] - norm**2 * aveQ[i]

stdMat = np.zeros(len(q)-1)
#for file in resultList:


t=np.arange(0,len(q)-1)
gamma=1.7

print("average position", norm*aveQ)

def theoDiff(x,a,b):
    return a*np.power(x,gamma)#+b

#startIndex = int(math.floor,t1/dt*0.5)
startIndex = int(np.floor(len(t)*0.5))
endIndex =int(len(t)-1)
#endIndex = int(np.floor(len(t)*0.24))
popt, pcov = curve_fit(theoDiff,t[startIndex:endIndex], varQ[startIndex:endIndex])
print(popt)





trajectory = plt.figure(1)
#plt.plot(trajectories)
plt.plot(q)
plt.xlabel('t')
plt.ylabel('sample trajectory')
trajectory.savefig("./plots/img/trajectory.pdf")

vQ = plt.figure(2)
plt.plot(t,varQ)
plt.plot(t[startIndex:endIndex],theoDiff(t[startIndex:endIndex],popt[0],popt[1]), color='#0066FF',linestyle='--',label=r'$\propto t^{\gamma}$')
plt.xlabel('t')
plt.ylabel('var(Q)')
vQ.savefig("./plots/img/varQ.pdf")
plt.legend()

vQlog = plt.figure(3)
plt.plot(t,varQ)
plt.plot(t[startIndex:endIndex],theoDiff(t[startIndex:endIndex],popt[0],popt[1]), color='#0066FF',linestyle='--',label=r'$\propto t^{\gamma}$')
plt.xscale('log', nonposx='clip')
plt.yscale('log', nonposy='clip')
plt.xlabel('t')
plt.ylabel('var(Q)')
vQlog.savefig("./plots/img/varQlog.pdf")




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
