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
    norm = np.power(len(resultList),-1.0)
    stdMat = np.zeros((len(q),len(resultList)))
    stdMat[:,0] = norm* norm* np.power(q[:,0],2) - norm*q[:,0]
    varQ = np.zeros(len(q))
    squaredQ = np.zeros(len(q))
    aveQ = np.zeros(len(q))

    for file in resultList:
        df = pd.read_csv(file)
        q = df.to_numpy()
        stdMat[:,resultList.index(file)] = norm*norm* np.power(q[:,0],2) - norm*q[0,:]
        squaQ = norm**2 * np.power(q,2)
        for i in range(0,len(squaredQ)):
            squaredQ[i] +=norm**2 * q[i]**2
            aveQ[i] += norm*q[i]

trajectories = np.zeros((len(q),6))
for i in range(0,6):
        random_index = randrange(len(resultList))
        df = pd.read_csv(resultList[i])
        m = df.to_numpy()
        trajectories[:,i]= norm*m[i,0]
   

stdErr = np.zeros(len(q))
for i in range(len(q)):
    stdErr[i] = np.var(stdMat[i,:])
print(stdErr)
 


for i in range(len(varQ)):
    varQ[i] =  squaredQ[i] -  aveQ[i]
varQ[0]=0

#tempArr = np.zeros(len(resultList))
#stdMat = np.zeros(len(q))
#for i in range(0,len(q)-1):
#    for file in resultList:
#        df = pd.read_csv(file)
#        q = df.to_numpy()
#        print(resultList.index(file))
#        tempArr[resultList.index(file)]=q[i]
#        stdMat[i]=np.var(tempArr)
        
        

t=np.arange(0,len(q))
gamma=1.9


def theoDiff(x,a,b):
    return a*np.power(x,b)

#startIndex = int(math.floor,t1/dt*0.5)
startIndex = int(np.floor(len(t)*0.4))
endIndex =int(len(t)-1)
#endIndex = int(np.floor(len(t)*0.43))
popt, pcov = curve_fit(theoDiff,t[startIndex:endIndex], varQ[startIndex:endIndex])
print(popt)

trajectory = plt.figure(1)
plt.plot(trajectories)
plt.plot(q)
plt.xlabel('t')
plt.ylabel('sample trajectory')
trajectory.savefig("./plots/img/trajectory.pdf")

errBarValue=50
skipValue =int(len(q)/errBarValue)

sli=slice(0,-1,skipValue)
tErr=t[sli]
varQErr=varQ[sli]
yerr=stdErr[sli]


logPlotStartIndex=int(0.01*len(q))
logsli=np.logspace(np.log10(logPlotStartIndex),np.log10(len(q)),errBarValue).astype(int)
tErrLog=t[logsli]
varQErrLog=varQ[logsli]
yerrLog=stdErr[logsli]


vQ = plt.figure(2)
#plt.errorbar(tErr,varQErr,yerr=yerr, ecolor='#FC9169')
plt.errorbar(tErr,varQErr,yerr=yerr,fmt='none', capsize=1.0, ecolor='#3B4CBF')
plt.plot(t,varQ,',', label='numerical results', color='#3B4CBF')
plt.plot(t[startIndex:endIndex],theoDiff(t[startIndex:endIndex],popt[0],popt[1]), color='#0066FF',linestyle='--',label=r'$\propto t^{\gamma}$')
plt.xlabel('t')
plt.ylabel('var(Q)')
plt.legend()
vQ.savefig("./plots/img/varQ.pdf")

vQlog = plt.figure(3)
plt.plot(t[logPlotStartIndex:-1],varQ[logPlotStartIndex:-1],',', label='numerical results', color='#3B4CBF')
#plt.errorbar(tErr,varQErr,yerr=yerr, ecolor='#FC9169')
plt.errorbar(tErrLog,varQErrLog,yerr=yerrLog,fmt='none', capsize=1.0, ecolor='#3B4CBF')
plt.plot(t[startIndex:endIndex],theoDiff(t[startIndex:endIndex],popt[0],popt[1]), linestyle='--',label=r'$\propto t^{\gamma}$')
#plt.xscale('log', nonposx='clip')
#plt.yscale('log', nonposy='clip')
plt.xscale('log' )
plt.yscale('log')


plt.xlabel('t')
plt.ylabel('var(Q)')
plt.legend()
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
