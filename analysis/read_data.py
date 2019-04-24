#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 13:54:46 2019

@author: joachim
"""

def read(filename, force_small = False):
    with open(filename, 'r') as file:
        number_of_bursts, p, full = [float(i.strip()) for i in file.readline().split(' ')]
        number_of_bursts, full = int(number_of_bursts), int(full)
        file_dump = file.readlines() #reads rest of the file..
        if force_small == False:
            if full == 1:
                bursts = [None]*number_of_bursts
                for line in range(number_of_bursts):
                    line_data = file_dump[line]
                    bursts[line] = [tuple([int(y) for y in j]) for j in [i[1:-1].split(" ") for i in line_data.split(', ')][:-1]]
        burst_sizes = file_dump[-1]
    
    burst_sizes = [int(i) for i in burst_sizes.split('\t')[0:number_of_bursts]]
    if force_small == False:    
        if full == 1:
            return burst_sizes, bursts
    return burst_sizes, []
    