#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 14:36:00 2019

@author: joachim
"""
import numpy as np
from read_data import read
import sys

def gamma_anal(sizes, multiplicity):
    if len(sizes)!=len(multiplicity):
        print('lenght doesnt match up...')
    avg_burst = sum([i*j for i,j in zip(sizes,multiplicity)])/sum(multiplicity)
    return avg_burst

def print_data(filename, x, y, header=None):
    if len(x) != len(y):
        print(' ...error! printing file {filename} with x,y not equal length!'.format(filename=filename))
    with open(filename, 'w') as file:
        if header:
            file.write(header)
        for i,j in zip(x,y):
            file.write('{x}\t{y}\n'.format(x=i, y=j))

def analyze_file(p, n, K, gamma_dict):
    freq_dict = {}
    
    for i in range(1,K+1):
        burst_sizes, bonds = read('../data/p{p}n{n}/run{i}p{p}n{n}.data'.format(i=i, p=p, n=n), force_small = True)
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
        N_dict = gamma_dict[N]
        p = np.array(list(N_dict.keys()))
        mean_clusters = np.array(list(N_dict.values())) #convert from python dict_type to np arrays
        print_data('./data/n{n}gamma'.format(n=N), p, mean_clusters)


P = [0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49]
N = [100000, 10000, 1000, 100, 10]
K = [1000]+[10000]*4

gamma_dict = {n:{} for n in N}
for n,k in zip(N,K):
    for p in P:        
        print('\nanalyzing {k} files with {n} bursts for probability {p}'.format(n=n, p=p, k=k))
        analyze_file(p, n, k, gamma_dict)
        
analyze_gamma_dict(gamma_dict)
