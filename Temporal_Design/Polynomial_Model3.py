# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 20:36:01 2022

@author: joeyv

[CODE EXPLANATION]
"""





"""
LIBRARIES
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


"""
CLASSES
"""
# Get random inputs for the independent variables
class getInput:
    
    def __init__(self):
        pass
    
    # Assigns a uniform random value to the provided range
    def getUniform(self, span):
        return np.random.uniform(span[0],span[1])
       
    # Assigns a normal random value based on mean and standard deviation
    def getNormal(self, mean, std_dev):
        return
    
# [Comment]
class getOutput1:
    
    # Organize inputs into self for the sequential solution
    def __init__(self):
        pass
        
    # Compute outputs for analysis 1
    def Analysis1(self, Path_vals, Path_vals_new):
        outs = np.array([2])
        ins = np.array([0, 1])
        Path_vals_new[2] = Path_vals[0] + Path_vals[1]
        return outs, ins
    
    # Compute outputs for analysis 2
    def Analysis2(self, Path_vals, Path_vals_new):
        outs = np.array([4, 5])
        ins = np.array([1, 2, 3])
        Path_vals_new[4] = Path_vals[2] - Path_vals[3]**2
        Path_vals_new[5] = 2*Path_vals[1] + Path_vals[3]
        return outs, ins
    
    # Compute outputs for analysis 3
    def Analysis3(self, Path_vals, Path_vals_new):
        outs = np.array([7])
        ins = np.array([4, 5, 6])
        Path_vals_new[7] = Path_vals[4] + math.sqrt(Path_vals[5]) - Path_vals[6]
        return outs, ins
    
    # Compute outputs for analysis 4
    def Analysis4(self, Path_vals, Path_vals_new):
        outs = np.array([9])
        ins = np.array([0, 7, 8])
        Path_vals_new[9] = Path_vals[7] - 2*(Path_vals[8] + Path_vals[0]**3)
        return outs, ins
    






"""
FUNCTIONS
"""





"""
USER INPUTS
"""
# Assign number of runs for each path
runs = 100

# Assign variable bounds
bounds = np.array([[1.0, 5.0], # x1
                   [1.0, 5.0], # x2
                   [0.0, 8.5], # x3
                   [0.5, 6.0], # x4
                   [-2.0, 9.0], # x5
                   [1.5, 10.0], # x6
                   [-3.5, 4.0], # x7
                   [-3.5, 2.5], # x8
                   [-4.0, 3.0], # x9
                   [1.0, 10.0]]) # x10

# List x-variable indices that will be independent
independ = np.array([[0, 1, 3, 6, 8], # Path1
                     [1, 2, 3, 6, 8], # Path2
                     [0, 4, 5, 6, 8], # Path3
                     [0, 1, 3, 7, 8]]) # Path4

# Define solver sequences
solver = getOutput1()
sequence = [[solver.Analysis1, solver.Analysis2,
             solver.Analysis3, solver.Analysis4], # Path1
            [solver.Analysis2, solver.Analysis1,
             solver.Analysis3, solver.Analysis4], # Path2
            [solver.Analysis3, solver.Analysis4,
             solver.Analysis1, solver.Analysis2], # Path3
            [solver.Analysis4, solver.Analysis1,
             solver.Analysis2, solver.Analysis3]] # Path4


"""
SCRIPT
"""

# Define matrix for path variables
Path_vals = np.zeros((np.shape(independ)[0],runs,np.shape(bounds)[0]))

# Define a vector for path variables that will be used to check inputs
Path_vals_new = np.zeros(np.shape(Path_vals)[2])

# Make class for creating random input
random = getInput()

# Update each run of each path with normalized random input within independent
## variable bounds
for i in range(0,np.shape(Path_vals)[0]):
    for j in range(0,runs):
        for k in range(0,np.shape(independ)[1]):
            index = independ[i,k]
            span = bounds[index,:]
            Path_vals[i,j,index] = random.getUniform(span)

# Define vectors for independent variable comparison
temp1 = np.zeros(np.shape(independ)[1])
temp2 = np.zeros(np.shape(independ)[1])

# [Solver - Comment]
for i in range(0,np.shape(Path_vals)[0]):
    for j in range(0,runs): 
        for k in range(0,np.shape(sequence)[1]):
                
            ## and maybe a while loop for the variable values to match up, with independent
            outs, ins = sequence[i][k](Path_vals[i,j,:],Path_vals_new)
            
        # Loop through the same run if old and new independ path vals do not match up
        for k in range(0,np.shape(independ)[1]):
            temp1[k] = Path_vals[i,j,independ[i,k]]
            temp2[k] = Path_vals_new[independ[i,k]]
        
            
            
        
        
        
        
# [Checker loop - Comment]



    




