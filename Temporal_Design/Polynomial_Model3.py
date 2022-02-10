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
    
# Solve for dependent outputs based on provided inputs
class getOutput:
    
    # Organize independent variables from dependent variables being solved for
    def __init__(self, independ, Path_vals):
        # Loop through a vector the size of independ assigning the Path_val of
        ## independent variables to correct x indexes and symbols for 
        ## dependent variables that are to be solved for.  Call it self.x to
        ## be used by whatever system is desired.  So x should have characters
        ## and floats...i.e. multiple data types
        pass
    
    # Solver for first system of equations
    def getSystem1(self):
        return
    
    # Solver for second system of equations
    def getSystem2(self):
        return
    
    
# Define nonlinear system of equations
#def func(x):
    #return [x[0] + x[1] - x[2],
            #x[2] - x[3]**2 - x[4],
            #2*x[1] + x[3] - x[5],
            #x[4] + math.sqrt(x[5]) - x[6] - x[7],
            #x[7] - 2*(x[8]+x[0]**3) - x[9]]




"""
FUNCTIONS
"""









"""
USER INPUTS
"""
# Assign number of runs for each path
runs = 10000

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





"""
SCRIPT
"""
# Define matrix for path variables
Path_vals = np.zeros((np.shape(independ)[0],runs,np.shape(bounds)[0]))

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
            
# Solve for dependent variable values


    




