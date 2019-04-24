#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 14:36:00 2019

@author: joachim
"""
import numpy as np
from read_file import read
import sys

def gamma_anal(sizes, multiplicity):
    if len(sizes)!=len(multiplicity):
        print('lenght doesnt match up...')
    avg_burst = sum([i*j for i,j in zip(sizes,multiplicity)])/len(sizes)
    return avg_burst
    

def analyze_file(p, n, K, gamma_dict):
    freq_dict = {}
    
    for i in range(1,K+1):
        burst_sizes, bonds = read('../data/p{p}n{n}/run{i}p{p}n{n}.data'.format(i=i, p=p, n=n), force_small = False)
        for size in burst_sizes:
            if size in freq_dict:
                freq_dict[size]+=1
            else:
                freq_dict[size]=1
        sys.stdout.write('\r ... completed file {i}'.format(i=i))

    burst_size = list(freq_dict.keys())
    multiplicity = list(freq_dict.values()) 
    mean_cluster = gamma_anal(burst_size, multiplicity)
    gamma_dict[n][p] = mean_cluster #this works..

def analyze_gamma_dict(gamma_dict):
    for N in gamma_dict: #Loop over N
        p = np.array(list(N.keys()))
        mean_clusters = np.array(list(N.values())) #convert from python dict_type to np arrays
        pc = 0.5 #for now.. might have to shift it to "remove finite lattice effects" (but lattice is effectively infinite..)
        plt.plot(np.log(p-pc), np.log(mean_clusters)) #this should be straight with exponent gamma. 

P = [0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49]
N = [100000]
K = [100]

gamma_dict = {n:{} for n in N}