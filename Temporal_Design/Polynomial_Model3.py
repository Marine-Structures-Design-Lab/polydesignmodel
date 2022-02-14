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
    def __init__(self,Path_vals):
        self.Path_vals = Path_vals
        pass
        
    # Compute outputs for analysis 1
    def Analysis1(self):
        self.Path_vals[2] = self.Path_vals[0] + self.Path_vals[1]
        return self.Path_vals
    
    # Compute outputs for analysis 2
    def Analysis2(self):
        self.Path_vals[4] = self.Path_vals[2] - self.Path_vals[3]**2
        self.Path_vals[5] = 2*self.Path_vals[1] + self.Path_vals[3]
        return self.Path_vals
    
    # Compute outputs for analysis 3
    def Analysis3(self):
        self.Path_vals[7] = self.Path_vals[4] + math.sqrt(self.Path_vals[5]) - self.Path_vals[6]
        return self.Path_vals
    
    # Compute outputs for analysis 4
    def Analysis4(self):
        self.Path_vals[9] = self.Path_vals[7] - 2*(self.Path_vals[8] + self.Path_vals[0]**3)
        return self.Path_vals
    




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

# [Solver - Comment]
for i in range(0,np.shape(Path_vals)[0]):
    for j in range(0,runs):
        
        # Make class for solver and establish path orders
        solver = getOutput1(Path_vals[i,j,:])
        sequence = [[solver.Analysis1, solver.Analysis2, solver.Analysis3, solver.Analysis4], # Path1
                    [solver.Analysis2, solver.Analysis1, solver.Analysis3, solver.Analysis4], # Path2
                    [solver.Analysis3, solver.Analysis4, solver.Analysis1, solver.Analysis2], # Path3
                    [solver.Analysis4, solver.Analysis1, solver.Analysis2, solver.Analysis3]] # Path4
        
        for k in range(0,np.shape(sequence)[1]):
                
            ## and maybe a while loop for the variable values to match up, with independent
            Path_vals[i,j,:] = sequence[i][k]()
            
            
            
        
        
        
        
# [Checker loop - Comment]



    




