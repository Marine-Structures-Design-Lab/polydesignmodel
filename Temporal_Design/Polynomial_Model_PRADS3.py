# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:36:50 2022

@author: joeyv

The goal of this code is to solve a sequential system of analyses for a
user-selected path.  The user is able to select the number of analyses, the
variables in the analyses, the equations in the analyses, the bounds on the
variables, among other things.  One caveat is that each analysis is only able
to be solved in one direction: inputs to outputs.  So iterations will be
required for situations where the inputs of one analysis are determined first,
when all or some of those inputs are also the outputs of an anlysis solved
later.  The user is able to select the type of minimizer, although BFGS or
CG methods are recommended.

On the first pass, all inputs are assigned a uniform random variable within
its accepted bounds.  Once the user is able to solve for all of the variables,
large rework loops take the most successful runs and reassign the inputs normal
random variables around those successful runs.  The simulation then runs again
for these new inputs.

The results focus on gathering the percentage of runs that each variable falls
within its required bounds, the percentage of runs having an alotted number of
successful variables, and the total number of small (L1) and large (L4) rework
loops.  Graphs are produced to display all of these results for each path.
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
            ### skip assignment if x is dependent var of analysis
            if (self.x[ind] in self.d):
                continue
            ### skip assignment if x has already been assigned as input
            elif (self.Pvn[ind] != 0):
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
                    expr[j] = expr[j].subs(self.x[ind],self.Pvn[ind])
        
        # Return the solvable function containing only one unknown variable
        return expr


# Assign the dependent output(s) solved for back to the numerical path vector
class assignOutput:
    
    # Initialize the class
    def __init__(self,sols,Path_vals_new,x):
        self.sols = sols
        self.Pvn = Path_vals_new;
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
            self.Pvn[ind_num] = self.sols[ind_var]
        
        # Return our updated numerical path vector
        return self.Pvn


# Determine if there any conflicts with the calculated outputs and inputs
class checkConflict:
    
    # Initialize the class
    def __init__(self,tol,Path_vals,Path_vals_new,depend,x):
        self.tol = tol
        self.Pv = Path_vals
        self.Pvn = Path_vals_new
        self.d = depend
        self.x = x
        return
    
    # Assess if calculated output value is adequately close to its input value
    def getCheck(self):
        
        # Create a list to store the True or False values
        checker = np.array(self.d,dtype=bool)
        
        # Loop through each solution of the analysis
        for i in range(0,np.shape(self.d)[0]):
            
            # Retrieve x variable index of the dependent variable
            ind = self.x.index(self.d[i])
            
            # Fill checker with the proper boolean value
            ### dependent variable not assigned an input
            if (self.Pv[ind] == 0):
                checker[i] = True
            ### dependent variable assignment and calculation are close
            elif (math.isclose(self.Pv[ind],self.Pvn[ind],rel_tol=self.tol)):
                checker[i] = True
            ### dependent variable assignment and calculation are not close
            else:
                checker[i] = False
        
        # return our boolean value(s)
        return checker


# Iterate through lower level loops if there are dependent variable conflicts
class getLoopy:
    
    # Initialize the class
    def __init__(self, Path_vals, Path_vals_new, Vars, analysis, depend, x, tol, l1_max, mini):
        self.Pv = Path_vals
        self.Pvn = Path_vals_new
        self.V = Vars
        self.analysis = analysis
        self.d = depend
        self.x = x
        self.tol = tol
        self.l1 = l1_max
        self.m = mini
        return
    
    # Do a rework loop w/in the analysis to try to get outputs to match inputs
    def analysisLoop(self):
        
        # Convert set to list
        self.V = list(self.V)
        
        # Initialize a list for the iterating x variable(s)
        x_iter = [];
        
        # Retrieve indices of variable in the analysis not assigned a value yet
        for i in range(0,np.shape(self.V)[0]):
            
            # Retrieve index for x variable in the analysis
            ind = self.x.index(self.V[i])
            
            # Determine if index value is unassigned
            if self.Pv[ind] == 0:
                x_list = [self.x[ind]]
                x_iter = x_iter + x_list
                
        # Initialize the check variable as false
        check = [False]*(np.shape(self.analysis)[0]);
            
        # Iterate over all equations of the analysis
        for i in range(0,np.shape(self.analysis)[0]):
            
            # Retrieve variables of the analysis's equation
            Vars = self.analysis[i].free_symbols
            
            # Iterate over all available x-iteration variables
            for j in range(0,np.shape(x_iter)[0]):
                
                # Check if equation has at least one iterating variable
                if x_iter[j] in Vars:
                    check[i] = True
                    break
            
        # Initialize L2-norm function
        l2norm = 0
        
        # Create symbolized L2-norm equation
        for i in range(0,np.shape(self.d)[0]):
            
            # Retrieve x variable index of the dependent variable
            ind = self.x.index(self.d[i])
            
            # Add term to function for dependent variable
            l2norm = l2norm + (self.analysis[i]/self.x[ind])**2
            
        # Add square root to the function
        l2norm = sp.sqrt(l2norm)
        
        # Create an L1 copy of the L2 norm expression
        l2norm_L1 = l2norm
        
        # If each equation has at least one iterating variable, perform L1 loop
        if all(check):
            
            # Initialize empty lists for optimization
            x0 = []
            xind = []
            Pvnind = []
            
            # Loop through variables in the set
            for i in range(0,np.shape(self.V)[0]):
                
                # Retrieve index for x variable in the analysis
                ind = self.x.index(self.V[i])
                
                # Substitute x variables in unless its an iterating variable
                if self.V[i] in x_iter:
                    x0 = x0 + [self.Pvn[ind]]
                    x_list = [self.x[ind]]
                    xind = xind + x_list
                    Pvnind = Pvnind + [ind]
                else:
                    l2norm_L1 = l2norm_L1.subs(self.x[ind],self.Pv[ind])
                
            # Optimize the L2-norm on the L1 rework loop
            l2norm_L1 = sp.utilities.lambdify(xind,l2norm_L1)
            ans_L1 = minimize(l2norm_L1,x0,method=self.m,tol=self.tol,options={'maxiter':self.l1})
            
            # Assign the new x value(s) to the new path values vector
            for i in range(0,np.shape(Pvnind)[0]):
                self.Pvn[Pvnind[i]] = ans_L1.x[i]
            
            # Retrieve the number of L1 loop iterations for the Rework matrix
            num_iters = ans_L1.nit
        
        else:
            
            # Set number of L1 iterations equal to 0
            num_iters = 0
            
        # Return new path values vector and number of iterations
        return self.Pvn, num_iters
    

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
            [2, 1, 3, 4], # Path2 and maybe 2,3,1,4 for comparison?
            [3, 4, 1, 2], # Path3 and variations?
            [4, 1, 2, 3]] # Path4

# Set tolerance for closeness of variables
tol = 1e-4

# Set max number of analysis (L1) rework loops
l1_max = 10

# Set max number of large (L4) rework loops
l4_max = 10

# Set method for the norm minimizer of the L1 rework loops
mini = 'BFGS'


"""
SCRIPT
"""
# Set up empty vectors and matrics
Path_vals = np.zeros((np.shape(sequence)[0],runs,np.shape(bounds)[0]))
Rework_L1 = np.zeros((np.shape(sequence)[0],runs,np.shape(sequence)[1]))
Rework_L4 = np.zeros((np.shape(sequence)[0],runs,1))

# Assign input values, calculate outputs, check for conflicts, resolve
for i in range(0,np.shape(Path_vals)[0]):        # i loops with paths
    for j in range(0,runs):                      # j loops with runs
        
        # Define a zero vector of path variables
        Path_vals_new = np.zeros(np.shape(Path_vals)[2])
        
        ###################### FIRST ANALYSIS #################################
        
        # Retrieve variables involved in first analysis
        index = sequence[i][0]
        variables = getVariables(analysis[index-1])
        Vars = variables.getVars()
        
        # Get random values for inputs of first analysis
        random = getInput(Vars,Path_vals_new,bounds,depend[index-1][:],x)
        Path_vals_new = random.getUniform()
        
        # Create function(s) for first analysis with numerical inputs and variable output(s)
        func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
        expr = func.getFunc()
        
        # Evaluate first analysis
        sols = sp.solve(expr)
        
        # Assign dependent variables(s) of first analysis to new path values vector
        solution = assignOutput(sols,Path_vals_new,x)
        Path_vals_new = solution.solAssign()
        
        # Assign a copy of new path values vector to official path values vector
        Path_vals[i,j,:] = np.copy(Path_vals_new)
        print(Path_vals[i,j,:])
        
        ######################### SECOND ANALYSIS #############################
        
        # Retrieve variables involved in the second analysis
        index = sequence[i][1]
        variables = getVariables(analysis[index-1])
        Vars = variables.getVars()
        
        # Get random values for inputs of second analysis
        random = getInput(Vars,Path_vals_new,bounds,depend[index-1][:],x)
        Path_vals_new = random.getUniform()
        
        # Create function(s) for second analysis with numerical inputs and variable output(s)
        func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
        expr = func.getFunc()
        
        # Evaluate second analysis
        sols = sp.solve(expr)
        
        # Assign dependent variables(s) of second analysis to new path values vector
        solution = assignOutput(sols,Path_vals_new,x)
        Path_vals_new = solution.solAssign()
        
        # Check for conflicts
        check = checkConflict(tol,Path_vals[i,j,:],Path_vals_new,depend[index-1][:],x)
        conflict = check.getCheck()
        
        # L1 loop if there are any conflicts
        if np.any(~conflict):
            
            # Gather new input values with the desired iterator
            # Populate L1 Rework loop with the number of iterations
            looper1 = getLoopy(Path_vals[i,j,:],Path_vals_new,Vars,analysis[index-1],depend[index-1][:],x,tol,l1_max,mini)
            Path_vals_new, Rework_L1[i,j,index-1] = looper1.analysisLoop()
            
            # Re-create function(s) for second analysis with numerical inputs and variable output(s)
            func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
            expr = func.getFunc()
            
            # Re-evaluate second analysis
            sols = sp.solve(expr)
            
            # Reassign dependent variable(s) of second analysis to new path values vector
            solution = assignOutput(sols,Path_vals_new,x)
            Path_vals_new = solution.solAssign()
            
            # THIS IS WHERE I NEED ANOTHER CONFLICT CHECKER AND TO DO A LARGE REWORK LOOP IF FAILED
            
        # Assign a copy of new path values vector to official path values vector
        Path_vals[i,j,:] = np.copy(Path_vals_new)
        print(Path_vals[i,j,:])
        
        ############################ THIRD ANALYSIS ###########################
        
        # Retrieve variables involved in the third analysis
        index = sequence[i][2]
        variables = getVariables(analysis[index-1])
        Vars = variables.getVars()
        
        # Get random values for inputs of third analysis
        random = getInput(Vars,Path_vals_new,bounds,depend[index-1][:],x)
        Path_vals_new = random.getUniform()
        
        # Create function(s) for third analysis with numerical inputs and variable output(s)
        func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
        expr = func.getFunc()
        
        # Evaluate third analysis
        sols = sp.solve(expr)
        
        # Assign dependent variables(s) of third analysis to new path values vector
        solution = assignOutput(sols,Path_vals_new,x)
        Path_vals_new = solution.solAssign()
        
        # Check for conflicts
        check = checkConflict(tol,Path_vals[i,j,:],Path_vals_new,depend[index-1][:],x)
        conflict = check.getCheck()
        
        # L1 loop if there are any conflicts
        if np.any(~conflict):
            
            # Gather new input values with the desired iterator
            # Populate L1 Rework loop with the number of iterations
            looper1 = getLoopy(Path_vals[i,j,:],Path_vals_new,Vars,analysis[index-1],depend[index-1][:],x,tol,l1_max,mini)
            Path_vals_new, Rework_L1[i,j,index-1] = looper1.analysisLoop()
            
            # Re-create function(s) for third analysis with numerical inputs and variable output(s)
            func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
            expr = func.getFunc()
            
            # Re-evaluate third analysis
            sols = sp.solve(expr)
            
            # Reassign dependent variable(s) of third analysis to new path values vector
            solution = assignOutput(sols,Path_vals_new,x)
            Path_vals_new = solution.solAssign()
            
            # THIS IS WHERE I NEED ANOTHER CONFLICT CHECKER AND TO DO A LARGE REWORK LOOP IF FAILED
            
        # Assign a copy of new path values vector to official path values vector
        Path_vals[i,j,:] = np.copy(Path_vals_new)
        print(Path_vals[i,j,:])
        
        ######################### FOURTH ANALYSIS #############################
        
        # Retrieve variables involved in the fourth analysis
        index = sequence[i][3]
        variables = getVariables(analysis[index-1])
        Vars = variables.getVars()
        
        # Get random values for inputs of fourth analysis
        random = getInput(Vars,Path_vals_new,bounds,depend[index-1][:],x)
        Path_vals_new = random.getUniform()
        
        # Create function(s) for fourth analysis with numerical inputs and variable output(s)
        func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
        expr = func.getFunc()
        
        # Evaluate fourth analysis
        sols = sp.solve(expr)
        
        # Assign dependent variables(s) of fourth analysis to new path values vector
        solution = assignOutput(sols,Path_vals_new,x)
        Path_vals_new = solution.solAssign()
        
        # Check for conflicts
        check = checkConflict(tol,Path_vals[i,j,:],Path_vals_new,depend[index-1][:],x)
        conflict = check.getCheck()
        
        # L1 loop if there are any conflicts
        if np.any(~conflict):
            
            # Gather new input values with the desired iterator
            # Populate L1 Rework loop with the number of iterations
            looper1 = getLoopy(Path_vals[i,j,:],Path_vals_new,Vars,analysis[index-1],depend[index-1][:],x,tol,l1_max,mini)
            Path_vals_new, Rework_L1[i,j,index-1] = looper1.analysisLoop()
            
            # Re-create function(s) for fourth analysis with numerical inputs and variable output(s)
            func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
            expr = func.getFunc()
            
            # Re-evaluate fourth analysis
            sols = sp.solve(expr)
            
            # Reassign dependent variable(s) of fourth analysis to new path values vector
            solution = assignOutput(sols,Path_vals_new,x)
            Path_vals_new = solution.solAssign()
            
            # THIS IS WHERE I NEED ANOTHER CONFLICT CHECKER AND TO DO A LARGE REWORK LOOP IF FAILED
            
        # Assign a copy of new path values vector to official path values vector
        Path_vals[i,j,:] = np.copy(Path_vals_new)
        print(Path_vals[i,j,:])

