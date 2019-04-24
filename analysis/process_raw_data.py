#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 15:16:25 2019

@author: joachim
"""

import sys
import numpy as np
from read_data import read

def plot_bins(bins):
    return [(bins[i]+bins[i+1])*(1/2) for i in range(len(bins)-1)]

def bin_middles(bins):
    return [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]

def exp_bins(max_num, base=2, nums=False):
    first_bin_edge = 50
    min_exp = int(np.log(first_bin_edge)/np.log(base))
    max_exp = int(np.ceil(np.log(max_num)/np.log(base)))
    bins = []
    if nums != False:
        max_exp = int(min_exp+nums)
    for i in range(min_exp, max_exp+1):
        bins.append(base**i)
    return bins

def filter_zeros(a,b):
    outa, outb = [], []
    for i,j in zip(a,b):
        if i != 0:
            outa.append(i)
            outb.append(j)
    return outa,outb
    
def add_frequencies(arr1,arr2): #returns the sum of arr1 and arr2 where the smallest is zero padded to be equal length to largest.
    arrays = [np.array(arr1), np.array(arr2)]
    index_biggest = int(len(arrays[0])<len(arrays[1])) #yeaheyahyeah not so readable... 
    index_smallest = int(not index_biggest)
    diff = len(arrays[index_biggest]) - len(arrays[index_smallest])
    return arrays[index_biggest] + np.hstack((arrays[index_smallest],np.zeros(shape=(diff,))))

def tau_b_anal(sizes, multiplicity):
    max_size = max(sizes)
    bin_spacing_log = np.exp(0.1) #keep spacing consistent...
    bins = exp_bins(max_size, base=bin_spacing_log)
    
    pdf, bins = np.histogram(sizes, bins = bins, weights = multiplicity)
    bins = plot_bins(bins)
    pdf, bins = filter_zeros(pdf, bins)
    
    cdf = np.cumsum(pdf)
    surv = np.max(cdf) - cdf
    surv, surv_bins = filter_zeros(surv, bins)    

    return (bins, pdf),(surv_bins, surv)

def fracdim_anal(bonds,spacing):
    flatten = [bond for burst in bonds for bond in burst]
    radius = [r(*bond[0:2]) for bond in flatten]
    bins = exp_bins(np.max(radius), base=spacing)

    f, bins = np.histogram(radius, bins=bins) #This is constant widht spacing in log space
    f,bins = filter_zeros(f,bins) #filter out zeros.
    
    return f, len(f)

def lin(x,a,b):
    return a*x+b

def pwr(x,a,b):
    return b*x**a

def r(x,y):
    return (x**2 + y**2)**0.5

def print_data(filename, x, y, header=None):
    if len(x) != len(y):
        print(' ...error! printing file {filename} with x,y not equal length!'.format(filename=filename))
    with open(filename, 'w') as file:
        if header:
            file.write(header)
        for i,j in zip(x,y):
            file.write('{x}\t{y}\n'.format(x=i, y=j))

def analyze_file(p, n, K):
    freq_dict = {}
    frac_dim_size_lim = 50
    frac_p_lim=0.465
    f_arr = np.zeros(shape=(10,))
    max_bin = 0

    for i in range(1,K+1):
        burst_sizes, bonds = read('../data/p{p}n{n}/run{i}p{p}n{n}.data'.format(i=i, p=p, n=n), force_small = False)
       
        if bonds and n>frac_dim_size_lim and p>frac_p_lim: #bonds not empty... (ant not n not (too small))
            freq, n_bins = fracdim_anal(bonds, frac_spacing) #If we dont need to check bond, we should definiatle force_small
            f_arr = add_frequencies(f_arr, freq)
            if n_bins > max_bin:
                max_bin = n_bins
                    
        for size in burst_sizes:
            if size in freq_dict:
                freq_dict[size]+=1
            else:
                freq_dict[size]=1
                
        sys.stdout.write('\r ... completed file {i}'.format(i=i))


    burst_size = list(freq_dict.keys())
    multiplicity = list(freq_dict.values())
    
    pdf, surv = tau_b_anal(burst_size, multiplicity)
    
    print_data('./data/p{p}n{n}{type}'.format(p=p, n=n, type='tau'), *pdf)
    print_data('./data/p{p}n{n}{type}'.format(p=p, n=n, type='survivor'), *surv)
    if bonds and n>frac_dim_size_lim and p>frac_p_lim:
        bins = plot_bins(exp_bins(1, base=frac_spacing, nums=max_bin))
        print_data('./data/p{p}n{n}{type}'.format(p=p, n=n, type='fracdim'), bins, f_arr)

P = [0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49]

N = [100000, 10000, 1000, 100, 10]
K = [1000]+[10000]*4

n_fits_s = []
n_fits_p = []

n_anal = len(P)*len(N)*K[0]
i_analyzed = 0
frac_spacing = np.exp(0.05)

for j,n in enumerate(N):
    k=K[j]
    for p in P:
        print('\nanalyzing {k} files with {n} bursts and probability {p}'.format(n=n, p=p, k=k))
        analyze_file(p, n, k)