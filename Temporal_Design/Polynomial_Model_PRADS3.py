# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:36:50 2022

@author: joeyv

[CODE EXPLANATION]
"""

"""
LIBRARIES
"""
import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve




"""
CLASSES
"""
# Get variables involved in the analysis equation(s)
class getVariables:
    
    # Initialize the class
    def __init__(self,analysis):
        self.analysis = analysis
        return
    
    # Find all variables involved in the analysis
    def getVars(self):
        
        # Initialize set
        Vars = {}
        
        # Retrieve variables of the analysis's equation(s)
        for i in range(0,np.shape(self.analysis)[0]):
            Vars[i] = self.analysis[i].free_symbols
        
        # Consolidate variables to one set without any repeats
        if np.shape(self.analysis)[0] > 1:
            for i in range(0,np.shape(self.analysis)[0]-1):
                Vars_new = Vars[i+1] - Vars[i]
                Vars[0] |= Vars_new
        
        # Return the joined set without any duplicates
        return Vars[0]


# Get random inputs for the independent variables
class getInput:
    
    # Initialize the class
    def __init__(self,Vars,Path_vals_new,bounds,depend,x):
        self.V = Vars
        self.Pvn = Path_vals_new
        self.bds = bounds
        self.d = depend
        self.x = x
        pass
    
    # Assign a uniform random value to the provided range
    def getUniform(self):
        
        # Convert set to list
        self.V = list(self.V)
        
        # Get uniform RV within bounds of particular x variable
        for i in range(0,np.shape(self.V)[0]):
            
            # Retrieve index for x variables in the analysis
            ind = self.x.index(self.V[i])
            
            # Get uniform random input within bounds if x is not dependent
            if self.x[ind] in self.d:
                continue
            else:
                self.Pvn[ind] = np.random.uniform(self.bds[ind,0],self.bds[ind,1])
        
        # Return vector with uniform RVs
        return self.Pvn
       
    # Assign a normal random value based on mean and standard deviation
    def getNormal(self, mean, std_dev):
        return


# Create an expression to be solved with numerical inputs replacing variable placeholders
class createFunction:
    
    # Initialize the class
    def __init__(self,analysis,Path_vals_new,Vars,depend,x):
        self.analysis = analysis
        self.Pvn = Path_vals_new
        self.V = Vars
        self.d = depend
        self.x = x
        return
    
    # Return a solvable function
    def getFunc(self):
        
        # Convert set to list
        self.V = list(self.V)
        print(self.V)
        
        # Make a copy of the analysis equations before any replacement
        expr = self.analysis
        
        # Loop through each variable of the equation
        for i in range(0,np.shape(self.V)[0]):
                
            # Retrieve index for x variables in the analysis
            ind = self.x.index(self.V[i])
            
            # Get uniform random input within bounds if x is not dependent
            if self.x[ind] in self.d:
                continue
            else:
                # Loop through each equation of the analysis
                for j in range(0,np.shape(self.analysis)[0]):
                    expr[j] = expr[j].subs(self.x[ind],self.Pvn[ind])
                    print(self.Pvn)
                    print(expr)
        
        # Return the solvable function containing only one unknown variable
        return expr


























# Evaluate function
class evalFuncs:
    
    # Initialize the class
    def __init__(self,Path_vals_new,analysis,depend,Vars,x):
        self.n = Path_vals_new
        self.analysis = analysis
        self.depend = depend
        self.Vars = Vars
        self.x = x
        return 
    
    # Return output of the desired analysis
    def evalAnalysis(self):
        
        # Convert set to list
        self.Vars = list(self.Vars)
        
        # Loop through each equation of the analysis
        for i in range(0,np.shape(self.analysis)[0]):
            
            # Loop through variables in an equation
            for j in range(0,np.shape(self.Vars)[0]):
                
                # Retrieve index for x variables in the analysis
                ind = self.x.index(self.Vars[j])
                
                # Create extra equations for independent variable values
                if (self.x[ind] in self.depend) or (i > 0):
                    continue
                else:
                    self.analysis.append(sp.Eq(self.x[ind]-self.n[ind]))
        print(self.n)
            
        # Solve the nonlinear system of equations
        print(self.Vars)
        answers = sp.solve(self.analysis, self.Vars)
        print(answers)
        
        
        
        
        
        
        
        
        
        
        # Add dependent x variable solution to vector
        #for i in range(0,np.shape(answers)[0]):
            #if (self.Vars[i] in self.depend):
                #self.n[i] = answers[i]
        
        #print(self.n)
        
        
        
        
        
        
        
        
        
        
        
        return self.n



"""
FUNCTIONS
"""





"""
USER INPUTS
"""
# Assign number of runs for each path
runs = 1

# Create symbols for all of the variables
x = sp.symbols('x1 x2 x3 x4 x5 x6 x7 x8 x9 x10')

# List out all of the equations
analysis = [[x[0] + x[1] - x[2]],                          # Analysis 1
            [x[2] - x[3]**2 - x[4], 2*x[1] + x[3] - x[5]], # Analysis 2
            [x[4] + sp.sqrt(x[5]) - x[6] - x[7]],          # Analysis 3
            [x[7] - 2*(x[8]+x[0]**3) - x[9]]];             # Analysis 4

# List out the variables being solved for in each analysis
depend = [[x[2]],       # Analysis 1
          [x[4], x[5]], # Analysis 2
          [x[7]],       # Analysis 3
          [x[9]]]       # Analysis 4

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

# Define solver sequences
sequence = [[1, 2, 3, 4], # Path1
            [2, 1, 3, 4], # Path2
            [3, 4, 1, 2], # Path3
            [4, 1, 2, 3]] # Path4












"""
SCRIPT
"""
# Set up empty vectors and matrics
Path_vals = np.zeros((np.shape(sequence)[0],runs,np.shape(bounds)[0]))
Rework = np.zeros((np.shape(sequence)[0],runs,np.shape(sequence)[1]))

# Assign input values, calculate outputs, check for conflicts, resolve
for i in range(0,np.shape(Path_vals)[0]):        # i loops with paths
    for j in range(0,runs):                      # j loops with runs
        
        # Define a zero vector of path variables
        Path_vals_new = np.zeros(np.shape(Path_vals)[2])
        
        # Retrieve variables involved in first analysis
        index = sequence[i][0]
        variables = getVariables(analysis[index-1])
        Vars = variables.getVars()
        #print(Vars)
        
        # Get random values for inputs of first analysis
        random = getInput(Vars,Path_vals_new,bounds,depend[index-1][:],x)
        Path_vals_new = random.getUniform()
        #print(Path_vals_new)
        
        # Create function(s) for first analysis with numerical inputs and variable output(s)
        func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
        expr = func.getFunc()
        print(expr)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        # Evaluate first analysis with fsolve
        #solver = evalFuncs(Path_vals_new,analysis[index-1][:],depend[index-1][:],Vars,x)
        #Path_vals_new = solver.evalAnalysis()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        # Get random numbers for inputs we do not have values for
        # [if Rework column 4 is 0...uniform input...else normal input]
        
        # Solve the equation for dependend variable(s)
        
        
        # Check -- Look for conflicts
        
        
        
        
        
        
        
        
        






