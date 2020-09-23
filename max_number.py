#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 10:28:11 2020

@author: andresnoe

OBDS Week 2 Day 2

Exercise 1 - Maximum number

Overall task: Find the maximum number
from an input list of numbers

1) Define list of numbers:
    26, 54, 93, 17, 77, 31, 44, 55, 20

2) Ask: Is element 1 greater than element 2?
    If yes: compare element 1 to element 3
        and repeat question
    If no: compare element 2 to element 3
        and repeat question
        
3) Repeat step until end 
"""


#Step 1
number_list = [26, 54, 93, 17, 77, 31, 44, 55, 20]

#Step 2.1: assume maximum number is initially
# the first number
number_max = number_list[0]

#Step 2.2: Work through for loop
for number in number_list:
    if number > number_max:
        number_max = number
#Step 3: Finish
print("The maximum number is", number_max)