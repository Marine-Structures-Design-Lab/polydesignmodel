# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 21:55:12 2022

@author: joeyv

[CODE EXPLANATION]
"""

"""
LIBRARIES
"""
import numpy as np
import math
import matplotlib.pyplot as plt



"""
CLASSES
"""
# Get random inputs for the independent variables
class getInput:
    
    # Initialize the class
    def __init__(self):
        pass
    
    # Assigns a uniform random value to the provided range
    def getUniform(self, span):
        return np.random.uniform(span[0],span[1])
       
    # Assigns a normal random value based on mean and standard deviation
    def getNormal(self, mean, std_dev):
        return


# Evaluate function
class evalFuncs:
    
    # Initialize the class
    def __init__(self, Path_vals):
        self.n = Path_vals
        return 
    
    # Return output of analysis 1
    def analysis1(self):
        self.n[2] = self.n[1] + self.n[2]
        return self.n
    
    # Return output of analysis 2
    def analysis2(self):
        self.n[4] = self.n[2] - self.n[3]**2
        self.n[5] = 2*self.n[1] + self.n[3]
        return self.n
    
    # Return output of analysis 3
    def analysis3(self):
        self.n[7] = self.n[4] + self.n[5]**(1/2) - self.n[6]
        return self.n
    
    # Return output of analysis 4
    def analysis4(self):
        self.n[9] = self.n[7] - 2*(self.n[8] + self.n[0]**3)
        return self.n


# Check if iterate and calculated value match within a tolerance
class checkIterate:
    
    # Initialize the class
    def __init__(self, Path_vals, iterate, func):
        self.Path_vals = Path_vals
        self.iterate = iterate
        self.func = func
        return
    
    # 
    
    
# Calculate a new iterate value for evaluation for desired nonlinear solver
class getIterate:
    
    # Initialize the class
    def __init__(self):
        pass
    
    # Newton's Method solver
    
    
    
    # Finite Difference method solver


"""
FUNCTIONS
"""
    









"""
USER INPUTS
"""
# Assign number of runs for each path
runs = 100

# Assign variable bounds
bounds = np.array([[1.0, 5.0],   # x1
                   [1.0, 5.0],   # x2
                   [0.0, 8.5],   # x3
                   [0.5, 6.0],   # x4
                   [-2.0, 9.0],  # x5
                   [1.5, 10.0],  # x6
                   [-3.5, 4.0],  # x7
                   [-3.5, 2.5],  # x8
                   [-4.0, 3.0],  # x9
                   [1.0, 10.0]]) # x10

# List x-variable indices that will be independent
independ = np.array([[0, 1, 3, 6, 8],  # Path1
                     [1, 2, 3, 6, 8],  # Path2
                     [0, 4, 5, 6, 8],  # Path3a
                     [0, 4, 5, 6, 8],  # Path3b
                     [0, 1, 3, 7, 8]]) # Path4

# Define solver sequences
sequence = [[1, 2, 3, 4], # Path1
            [2, 1, 3, 4], # Path2
            [3, 4, 1, 2], # Path3a
            [3, 4, 1, 2], # Path3b
            [4, 1, 2, 3]] # Path4
for i in range(np.shape(sequence)[0]):
    for j in range(np.shape(sequence)[1]):
        sequence[i][j] = 'solver.analysis' + str(sequence[i][j])
        
# List iterating variable(s) of each path



"""
SCRIPT
"""
# Find x-variable indices that will be dependent
depend = np.zeros((np.shape(independ)[0],np.shape(independ)[1]))
for i in range(0,np.shape(independ)[0]):
    count = 0
    for j in range(0,np.shape(bounds)[0]):
        if j in independ[i,:]:
            continue
        else:
            depend[i,count] = j
            count += 1
del count

# Set up empty vectors and matrics
Path_vals = np.zeros((np.shape(independ)[0],runs,np.shape(bounds)[0]))
Path_vals_new = np.zeros(np.shape(Path_vals)[2])
Rework = np.zeros((np.shape(independ)[0],runs,np.shape(sequence)[1]))

# Create random input class
random = getInput()

# Populate Path_vals with uniform RVs for independent variables of each path
for i in range(0,np.shape(Path_vals)[0]):        # i loops with paths
    for j in range(0,runs):                      # j loops with runs
        for k in range(0,np.shape(independ)[1]): # k loops with independent variables
            index = independ[i,k]
            span = bounds[index,:]
            Path_vals[i,j,index] = random.getUniform(span)

# Populate Path_vals with solutions for dependent variables of each path
for i in range(0,np.shape(Path_vals)[0]):        # i loops with paths
    for j in range(0,runs):                      # j loops with runs
        for k in range(0,np.shape(sequence)[1]): # k loops with analyses
            
            # Create analysis class and evaluate for new Path values
            solver = evalFuncs(Path_vals[i,j,:])
            Path_vals_new = eval(sequence[i][k] + '()')
            
            # while loop
            
            
            # Replace Path values with new Path values for a run
            Path_vals[i,j,:] = Path_vals_new
            
            




















