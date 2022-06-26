# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 16:57:40 2022

@author: joeyv

Create utility theory code for the set-space reduction of the polynomial model
and track the effectiveness of the space reductions for a set number of Monte
Carlo runs of each path (and all analyses collectively)
"""

"""
LIBRARIES
"""
import numpy as np
import sympy as sp
import math
from scipy.optimize import minimize
import matplotlib.pyplot as plt

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
    def __init__(self,Vars,Path_vals,depend,x):
        self.V = Vars
        self.Pv = Path_vals
        self.d = depend
        self.x = x
        return
    
    # Assign a uniform random value to the provided range
    def getUniform(self, sbounds):
        
        # Convert set to list
        self.V = list(self.V)
        
        # Loop through the x variables of the analysis
        for i in range(0,np.shape(self.V)[0]):
            
            # Retrieve index for x variable in the analysis
            ind = self.x.index(self.V[i])
            
            # Skip assignment if x is dependent variable of analysis
            if (self.x[ind] in self.d):
                continue
            # Skip assignment if x has already been assigned as input
            elif (self.Pv[ind] != 0):
                continue
            # Get uniform random input within bounds
            else:
                self.Pv[ind] = np.random.uniform(sbounds[ind,0],sbounds[ind,1])
        
        # Return vector with uniform RVs
        return self.Pv


# Create an expression to be solved with numerical inputs replacing variable placeholders
class createFunction:
    
    # Initialize the class
    def __init__(self,analysis,Path_vals,Vars,depend,x):
        self.analysis = analysis
        self.Pv = Path_vals
        self.V = Vars
        self.d = depend
        self.x = x
        return
    
    # Return a solvable function
    def getFunc(self):
        
        # Convert set to list
        self.V = list(self.V)
        
        # Make a copy of the analysis equations before any replacement
        expr = np.copy(self.analysis)
        
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
                    expr[j] = expr[j].subs(self.x[ind],self.Pv[ind])
        
        # Return the solvable function containing only one unknown variable
        return expr


# Assign the dependent output(s) solved for back to the numerical path vector
class assignOutput:
    
    # Initialize the class
    def __init__(self,sols,Path_vals,x):
        self.sols = sols
        self.Pv = Path_vals;
        self.x = x
        return
    
    # Assign the calculated solution(s) to the path vector
    def solAssign(self):
        
        # Loop through all solutions of the analysis
        for i in range(0,len(self.sols)):
        
            # Retrieve x variable index of the solution
            ind_var = list(self.sols.keys())[i]
            
            # Retrieve numerical index of the x variable index
            ind_num = self.x.index(ind_var)
            
            # Add solution to the proper index in the numerical path vector
            num = self.sols[ind_var]
            if sp.im(num) != 0:
                self.Pv[ind_num] = np.NaN
            else:
                self.Pv[ind_num] = self.sols[ind_var]
        
        # Return our updated numerical path vector
        return self.Pv




"""
FUNCTIONS
"""
# Check success of variables during simulation
def solChecker(Path,bounds,Vars,j,x):
    
    # Convert set to list
    Vars = list(Vars)
    
    # Gather success results for the analysis run
    for i in range(0,len(Vars)):
        
        # Retrieve index for x variable in the analysis
        ind = x.index(Vars[i])

        # Check if variable value is within its bounds
        if (Path[ind] >= bounds[ind,0]) and (Path[ind] < bounds[ind,1]):
            continue 
        else:
            return False
    
    # Return true for successful analysis variables
    return True




"""
USER INPUTS
"""
# Assign maximum number of runs for each path
runs = 10

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
            [2, 3, 1, 4], # Path2
            [3, 4, 1, 2], # Path3
            [4, 1, 2, 3]] # Path4








"""
SCRIPT
"""
# Track set-based design space
Set_bounds = np.zeros((1,np.shape(bounds)[0],np.shape(bounds)[1]))
Set_bounds[0,:,:] = np.copy(bounds)

# Create arrays for Monte Carlo run values
Path_vals = np.zeros((runs,len(x)))
Path_vals_s = np.zeros((1,len(x)))

# Set iteration counter
count = 0

# Loop through each analysis of path one once
for i in range(0,np.shape(sequence)[1]):

    # Append bounds from previous iteration
    Set_bounds = np.append(Set_bounds,np.copy(Set_bounds[count,:,:].reshape((1,np.shape(bounds)[0],np.shape(bounds)[1]))),axis=0)
    
    # Retrieve variables involved in analysis
    index = sequence[0][i]
    variables = getVariables(analysis[index-1])
    Vars = variables.getVars()
    
    # Loop through each Monte Carlo run
    for j in range(0,np.shape(Path_vals)[0]):
        
        # Get random values for inputs of analysis
        random = getInput(Vars,Path_vals[j,:],depend[index-1][:],x)
        Path_vals[j,:] = random.getUniform(Set_bounds[count,:,:])
        
        # Create function(s) for analysis with numerical inputs and variable output(s)
        func = createFunction(analysis[index-1],Path_vals[j,:],Vars,depend[index-1][:],x)
        expr = func.getFunc()
        
        # Evaluate analysis
        sols = sp.solve(expr)
        
        # Assign dependent variables(s) of analysis to new path values vector
        solution = assignOutput(sols,Path_vals[j,:],x)
        Path_vals[j,:] = solution.solAssign()
        
        # Check if all variables fall within bounds
        check = solChecker(Path_vals[j,:],bounds,Vars,j,x)
        
        # Collect successful analysis indices
        if check and np.sum(Path_vals_s)==0:
            Path_vals_s[0,:] = np.copy(Path_vals[j,:])
        elif check:
            Path_vals_s = np.append(Path_vals_s,[Path_vals[j,:]],axis=0)
    
    # Find min/max variable values of collected indices
    
    
    # [Produce histogram of viable design space for each variable - This could be used as experience later if desired instead of just cutting off zero utility values]
    
    # Change variable set bounds based on collected indices
    
    # Reset Path values to zeros
    Path_vals = np.zeros((runs,len(x)))
    Path_vals_s = np.zeros((1,len(x)))
        
    # Increase iteration counter by 1
    count = count + 1

# Plot upper and lower bounds for each analysis iteration of the path
fig = plt.figure(figsize=(10, 6))
xi = np.arange(np.shape(Set_bounds)[0])
color = iter(plt.cm.rainbow(np.linspace(0, 1, np.shape(bounds)[0])))
for i in range(0,np.shape(Set_bounds)[1]):
    c = next(color)
    for j in range(0,np.shape(Set_bounds)[2]):
        if j == 0:
            plt.plot(xi,Set_bounds[:,i,j],'-o',c=c,label=x[i])
        else:
            plt.plot(xi,Set_bounds[:,i,j],'-o',c=c)
plt.title('Path 1')
plt.xlabel('Analysis Iteration')
plt.ylabel('Bound')
plt.legend(loc='best')
plt.grid(which='major',axis='both')
plt.show()


# Create RVs for analysis inputs based on Set_bounds
# Calculate analysis output(s)






# All together

# Start with set space intervals based on the input bounds
# Create random values for inputs of analysis 1 within their bounds
# Use those values to calculate x3
# Plot PDF of values for x1 and x2 that allow x3 to pass
# Reduce design space for intervals having 0% of values
# Do same for other 3 analyses while uninfluenced by other info from others
# Return new set space intervals























