#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 13:34:33 2020

@author: andresnoe

Bubble sort a list of numbers

Pseudocode
for i from 1 to N:
    for j from 0 to N-1:
        if a[j]>a[j+1]
        swap(a[j],a[j+1])

"""

number_list = [26, 54, 93, 17, 77, 31, 44, 55, 20]

for j in range(len(number_list),1,-1): # STEPPING LOOP
    # Use the number j which starts at length of list and goes down to one in increments of -1
    # Repeat the SWAPPING LOOP the required numnber of times
    for i in range(0,j-1,1): #SWAPPING LOOP
    # Compare two numbers and swap 
        if number_list[i]>number_list[i+1]: # Pairwise comparisons
            # SWAP THE TWO NUMBERS IF TRUE
            temp_variable = number_list[i] # Create temporary variable
            number_list[i] = number_list[i+1] # Copy/paste one element to have a duplicate in new position and old position
            number_list[i+1] = temp_variable # Use temporary variable to write over old position
        print(number_list)

# Best practice: rename variables to have meaningful variable names
# Eg j = step/iteration, i = current index