#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 15:16:25 2019

@author: joachim
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from read_data import read

def plot_bins(bins):
    return [(bins[i]+bins[i+1])*(1/2) for i in range(len(bins)-1)]

def exp_bins(max_num, base=2):
    first_bin_edge = 50
    min_exp = int(np.log(first_bin_edge)/np.log(base))
    max_exp = int(np.ceil(np.log(max_num)/np.log(base)))
    bins = []
    for i in range(min_exp, max_exp+1):
        bins.append(base**i)
    return bins

def filter_zeros(a,b):
    outa, outb = [], []
    for i,j in zip(a,b):
        if i!=0:
            outa.append(i)
            outb.append(j)
    return outa,outb
    
def survivor_plot(sizes, multiplicity, plot=False):
    max_size = max(sizes)
    bin_spacing_log = np.exp(0.1) #keep spacing consiteng...
    bins = exp_bins(max_size, base=bin_spacing_log)
    
    pdf, bins = np.histogram(sizes, bins = bins, weights = multiplicity)
    bins = plot_bins(bins)
    pdf, bins = filter_zeros(pdf, bins)
    
    cdf = np.cumsum(pdf)
    surv = np.max(cdf) - cdf
    surv, surv_bins = filter_zeros(surv, bins)    

    fitrange=slice(10,50)
    
    survivor_fit, pdf_fit, bins_fit, surv_bins_fit = np.log(surv[fitrange]), np.log(pdf[fitrange]), np.log(bins[fitrange]), np.log(surv_bins[fitrange])


    pfit, pcov = opt.curve_fit(lin, bins_fit, pdf_fit)
    sfit, scov = opt.curve_fit(lin, surv_bins_fit, survivor_fit)
    
    if plot:
        plt.loglog(bins,surv)
        plt.loglog(bins,pdf)
        plt.loglog(bins,cdf)
            
        plt.plot(np.exp(bins_fit), [pwr(x, pfit[0], np.exp(pfit[1])) for x in np.exp(bins_fit)],'--')
        plt.plot(np.exp(surv_bins_fit), [pwr(x, sfit[0], np.exp(sfit[1])) for x in np.exp(bins_fit)],'--')
        
        plt.legend(['survivor function', 'PDF', 'CDF'])
        plt.title('Tau fit and survivor function fit.\n Tau={t:.04f}, Tau_survivor={ts:.04f}. b=3/2*Tau_survivor={b:.04f}'.format(t=-pfit[0], ts=-sfit[0], b=-sfit[0]*3/2))
        plt.xlabel('Ln(s)')
        plt.ylabel('Ln(N_s)')

    return -pfit[0], -sfit[0]

def lin(x,a,b):
    return a*x+b

def pwr(x,a,b):
    return b*x**a

P = [0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49]
N = [10000]
k = 250

n_fits_s = []
n_fits_p = []

n_anal = len(P)*len(N)*k
i_analyzed = 0

for n in N:
    fits_s = []
    fits_p = []
    for p in P:
        freq_dict = {}
        for i in range(1,k+1):
            i_analyzed += 1
            sys.stdout.write('\r {i:.01f}% Completed'.format(i=100*i_analyzed/n_anal))
            burst_sizes = read('../data/p{p}n{n}/run{i}p{p}n{n}.data'.format(i=i, p=p, n=n), force_small = True)
            for size in burst_sizes:
                if size in freq_dict:
                    freq_dict[size]+=1
                else:
                    freq_dict[size]=1
    
        burst_size = list(freq_dict.keys())
        multiplicity = list(freq_dict.values())
        pfit, sfit = survivor_plot(burst_size, multiplicity)
        fits_s.append(sfit)
        fits_p.append(pfit)
        plt.close('all')

    n_fits_s.append(fits_s)
    n_fits_p.append(fits_p)

for nfits,pfits,size in zip(n_fits_s,n_fits_p,N):
    plt.plot(P,nfits,'o', label='{} Bursts survivor'.format(size))
    plt.plot(P,pfits,'o', label='{} Burst tau'.format(size))
    
    fitfit_s,covs = opt.curve_fit(lin, P, nfits)
    fitfit_p,covp = opt.curve_fit(lin, P, pfits) #linear fits

    plt.plot(P, [lin(x,*fitfit_s) for x in P],'g')
    plt.plot(P, [lin(x,*fitfit_p) for x in P],'r') #in the data region
    
    p_extrap=[P[-1], 0.5]
    
    plt.plot(p_extrap, [lin(x,*fitfit_s) for x in p_extrap],'g--')
    plt.plot(p_extrap, [lin(x,*fitfit_p) for x in p_extrap],'r--') #extraploated.
    

    

plt.legend()