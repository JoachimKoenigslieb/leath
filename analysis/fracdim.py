#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:05:14 2019

@author: joachim
"""

import numpy as np
import matplotlib.pyplot as plt
from read_data import read
import scipy.optimize as opt

def r(x,y):
    return (x**2 + y**2)**0.5

def fitfunc(x, a, b):
    return x**a+b

def lin(x, a, b):
    return a*x+b

def bin_middles(bins):
    return [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]

def exp_bins(max_num, base=2):
    first_bin_edge = 50
    min_exp = int(np.log(first_bin_edge)/np.log(base))
    max_exp = int(np.ceil(np.log(max_num)/np.log(base)))
    bins = []
    for i in range(min_exp, max_exp+1):
        bins.append(base**i)
    return bins

fracdim = []

n_bins = 100

f_arr = np.zeros((n_bins,))
for i in range(1, 11):
    print(i)
    freq, bonds = read('../data/p0.49n10000/run'+str(i)+'p0.49n10000.data')
    flatten = [bond for burst in bonds for bond in burst]
    radius = [r(*bond[0:2]) for bond in flatten]
    f, bins = np.histogram(radius, bins=n_bins)
    f_arr += f    
bins = bin_middles(bins)

M = 1/10000*np.cumsum(f_arr) #get mass as the integral of f_arr. 
logm = np.log(M)
logr = np.log(bins) 

s=slice(0,-1)

logfit, logcov = opt.curve_fit(lin, logr[s], logm[s]) #filter end and start
plt.plot(logr, logm,'o')
plt.plot(logr[s], [lin(x, *logfit) for x in logr[s]], '--')
plt.title('10000 clusters of size 10000. Fractal dimension is d={:.04f}'.format(logfit[0]))
plt.xlabel('ln(r)')
plt.ylabel('avg. ln(M)')
plt.show()