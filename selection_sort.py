#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 10:57:48 2020

@author: andresnoe

OBDS Week 2 Day 2

Exercise 2 - Sort numbers

Overall task: Sort numbers in ascending order
from an input list of numbers

1) Define list of numbers:
    26, 54, 93, 17, 77, 31, 44, 55, 20

2) Ask: Is element 1 greater than element 2?
    If yes: compare element 1 to element 3
        and repeat question
    If no: compare element 2 to element 3
        and repeat question

"""

number_list = [26, 54, 93, 17, 77, 31, 44, 55, 20]
def max_number (number_list):
    number_max = number_list[0]
    max_index = 0
    loop_index = 0
    
    for number in number_list:
        if number > number_max:
            number_max = number
            max_index = loop_index
        # print(number, loop_index)
        loop_index += 1 

    return [number_max, max_index] 
# print(max_number(number_list))
    
# the below codes sort the number_list in ascending order

for i in range(len(number_list),1,-1):
    biggest_number, max_index = max_number(number_list[:i])
    temp_variable = number_list[i-1]
    number_list[i-1] = biggest_number
    number_list[max_index] = temp_variable
    print(number_list[:i])
    print (i, biggest_number, max_index)
    
print(number_list)