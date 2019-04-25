#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:23:40 2019

@author: joachim
"""

import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
import numpy as np

def getp(filename):
    n_index = filename.index('n')
    return float(filename[1:n_index])

def getn(filename):
    n_index = filename.index('n') #finds the first 'n'. Walks trough stirng untill finds a not number
    for j,c in enumerate(filename[n_index+1:]): #holy shit who writes this kind of shit code..?
        if c.isalpha():
            break
    return int(filename[n_index+1:n_index+j+1])

def pairsort(l1,l2):
    return list(zip(*sorted(zip(l1,l2), key = lambda item: item[0]))) #sort from the first index

def sort_n(d, files):
    num_set = set()
    for f in files: #lets get all N's in the dataset.
        num = getn(f)
        if num not in num_set:
            num_set.add(int(num))
            d[num] = []
        d[num].append(f)
    return d

def read(filename):
    with open('./data/{}'.format(filename), 'r') as file:
        x,y = [],[]
        for line in file:
            i,j = [float(k) for k in line.strip('\n').split('\t')]
            x.append(i)
            y.append(j)
    return x,y

def lin(x,a,b):
    return a*x+b

def pwr(x,a,b):
    return b*x**a

def linfit(x,y,s, name=None):
    fit, cov = curve_fit(lin, x[s], y[s])
    if name != None:
        plt.ioff()
        f,a = plt.subplots(1,1)
        a.plot(x,y,'ok')
        a.plot(x[s], [lin(i,*fit) for i in x[s]], 'r--') #fitregion 
        a.set_title('Slope: {s:.03f}'.format(s=fit[0]))
        f.savefig('./data/fits/{name}'.format(name=name))
        plt.close(f)
        plt.ion()
    return fit,cov 

def sort_legend(ax):
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: int(t[0].split(' ')[0])))
    return handles, labels

files = os.listdir('./data')

survivor, tau, fracdim, gamma = [],[],[],[]
for f in files:
    if 'tau' in f:
        tau.append(f)
    if 'survivor' in f:
        survivor.append(f)
    if 'fracdim' in f:
        fracdim.append(f)
    if 'gamma' in f:
        gamma.append(f)

tau_dict=sort_n({},tau)
surv_dict=sort_n({},survivor)
fracdim_dict=sort_n({},fracdim)
gamma_dict=sort_n({},gamma)

plt.ioff()

tau_fits = {}
surv_fits = {}
fracdim_fits = {}
gamma_fits = {}

plt.close('all')

markers = ['o','s','x','^','d']

save_fits = False

for N in tau_dict:
    filenames = tau_dict[N]
    P = [getp(f) for f in filenames] #extract P
    s=slice(2,25) #pretty arbitrary but we can use it for now..    
    if save_fits:
        tau = [-linfit(*[np.log(i) for i in read(f)],s, name=f+'.png')[0][0] for f in filenames] 
    else:
        tau = [-linfit(*[np.log(i) for i in read(f)],s, name=None)[0][0] for f in filenames] 
    tau_fits[N] = pairsort(P,tau)
    
for N in surv_dict:
    filenames = surv_dict[N]
    P = [getp(f) for f in filenames] 
    s=slice(2,25) #pretty arbitrary but we can use it for now..
    if save_fits:
        surv_tau = [-linfit(*[np.log(i) for i in read(f)],s, name=f+'.png')[0][0] for f in filenames] 
    else:
        surv_tau = [-linfit(*[np.log(i) for i in read(f)],s, name=None)[0][0] for f in filenames] 
    b=[3/2*i for i in surv_tau]
    surv_fits[N] = pairsort(P,b)

for N in fracdim_dict:
    filenames = fracdim_dict[N]
    P = [getp(f) for f in filenames]
    s=slice(0,15) #pretty arbitrary but we can use it for now..
    if save_fits:
        fracdim = [linfit(*[np.log(i) for i in read(f)],s, name=f+'.png')[0][0] for f in filenames] 
    else:
        fracdim = [linfit(*[np.log(i) for i in read(f)],s, name=None)[0][0] for f in filenames] 
    fracdim_fits[N] = pairsort(P,fracdim)

for N in gamma_dict:
    f = gamma_dict[N][0]
    P, mean_clusters = read(f)
    P = np.log(abs(0.5-np.array(P)))
    mean_clusters = np.log(mean_clusters)
    if save_fits:
        gamma = linfit(P, mean_clusters, slice(0,None), name=f+'.png')[0][0]
    else:
         gamma = linfit(P, mean_clusters, slice(0,None))[0][0]
    gamma_fits[N] = gamma

f_surv, a_surv = plt.subplots(1,1)
for N,m in zip(surv_fits,markers):
    P, surv = surv_fits[N]
    a_surv.plot(*surv_fits[N],'--k', marker=m, label='{n} bursts'.format(n=N))
#    a_surv.loglog(abs(np.array(P)-0.5), surv, '--k', marker=m, label='{n} bursts'.format(n=N))
#a_surv.set_title('Survival fit')
a_surv.set_xlabel('p')
a_surv.set_ylabel('b')
a_surv.legend(*sort_legend(a_surv))
    
f_tau, a_tau = plt.subplots(1,1)
for N,m in zip(tau_fits,markers):
    P, tau = tau_fits[N]
    a_tau.plot(*tau_fits[N],'--k', marker=m, label='{n} bursts'.format(n=N))
#    a_tau.loglog(abs(np.array(P)-0.5), tau, '--k', marker=m, label='{n} bursts'.format(n=N))
#a_tau.set_title('Tau fit')
a_tau.set_xlabel('p')
a_tau.set_ylabel('tau')
a_tau.legend(*sort_legend(a_tau))

f_frac, a_frac = plt.subplots(1,1)
for N,m in zip(fracdim_fits,markers): 
    if N>100:
        P, d = fracdim_fits[N]
        a_frac.plot(P,d, color='k',linestyle="None", marker=m, label='{n} bursts'.format(n=N))
        fracdim_extrap_fit = linfit(P,d,slice(0,None))[0]
        a_frac.plot(P, [lin(x, *fracdim_extrap_fit) for x in P],'k')
        P_extrap = [P[-1],0.5]
        a_frac.plot(P_extrap, [lin(x, *fracdim_extrap_fit) for x in P_extrap], 'k--')
        
#a_frac.set_title('Fractal dimension fits')
a_frac.set_xlabel('p')
a_frac.set_ylabel('d')
a_frac.legend()

f_gamma, a_gamma = plt.subplots(1,1)
a_gamma.semilogx(*pairsort(list(gamma_fits.keys()), list(gamma_fits.values())),'ok--', label='Gamma for different values of N')
#a_gamma.set_title('Gamma fits semilogx scale')
a_gamma.set_xlabel('Number of bursts')
a_gamma.set_ylabel('gamma')
a_gamma.legend()


plt.show()
