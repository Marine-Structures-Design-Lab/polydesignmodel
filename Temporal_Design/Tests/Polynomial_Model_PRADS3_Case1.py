# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:36:50 2022

@author: Joseph B. Van Houten (joeyvan@umich.edu)

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

User requirements: Should only need to make changes to the 'USER INPUTS' 
section and the graph-related code at the end of the 'SCRIPT' section.
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
    def __init__(self,Vars,Path_vals_new,depend,x):
        self.V = Vars
        self.Pvn = Path_vals_new
        self.d = depend
        self.x = x
        pass
    
    # Assign a uniform random value to the provided range
    def getUniform(self, bounds):
        
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
            elif (self.Pvn[ind] != 0):
                continue
            # Get uniform random input within bounds
            else:
                self.Pvn[ind] = np.random.uniform(bounds[ind,0],bounds[ind,1])
        
        # Return vector with uniform RVs
        return self.Pvn
       
    # Assign a normal random value based on mean and standard deviation
    def getNormal(self, mean, std, bounds, mv, stdr, h):
        
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
            elif (self.Pvn[ind] != 0):
                continue
            # Get normal random input given mean and standard deviation
            else:
                if (mv[ind] == 0) or (np.any(np.isnan(mv))):
                    self.Pvn[ind] = np.random.normal(mean[ind],std[ind])
                else:
                    std_dev = stdr*(np.abs(bounds[ind,1]-bounds[ind,0]))/(0.5**(h-1))
                    self.Pvn[ind] = np.random.normal(mean[ind],std_dev)
        
        # Return vector with normal RVs
        return self.Pvn
    
    # Assign a combination of uniform and normal random variables
    def getHybrid(self, bounds, mean, stdr, h):
        
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
            elif (self.Pvn[ind] != 0):
                continue
            # Get uniform or normal input
            else:
                # Assign uniform RV if we do not have a mean value
                if (mean[ind] == 0) or (np.any(np.isnan(mean))):
                    self.Pvn[ind] = np.random.uniform(bounds[ind,0],bounds[ind,1])
                # Assign normal RV with diminishing standard deviation if we do have a mean value
                else:
                    std_dev = stdr*(np.abs(bounds[ind,1]-bounds[ind,0]))/(0.5**(h-1))
                    self.Pvn[ind] = np.random.normal(mean[ind],std_dev)
        
        # Return vector with uniform and normal RVs
        return self.Pvn


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
            num = self.sols[ind_var]
            if sp.im(num) != 0:
                self.Pvn[ind_num] = np.NaN
            else:
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
            elif (math.isclose(self.Pvn[ind],self.Pv[ind],rel_tol=self.tol)):
                checker[i] = True
            ### dependent variable assignment and calculation are not close
            else:
                checker[i] = False
        
        # return our boolean value(s)
        return checker


# Iterate through lower level loops if there are dependent variable conflicts
class getLoopy:
    
    # Initialize the class
    def __init__(self, Path_vals, Path_vals_new, Vars, analysis, depend, x, l1_max, mini):
        self.Pv = Path_vals
        self.Pvn = Path_vals_new
        self.V = Vars
        self.analysis = analysis
        self.d = depend
        self.x = x
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
            ans_L1 = minimize(l2norm_L1,x0,method=self.m,options={'maxiter':self.l1})
            
            # Assign the new x value(s) to the new path values vector
            for i in range(0,np.shape(Pvnind)[0]):
                self.Pvn[Pvnind[i]] = ans_L1.x[i]
            
            # Retrieve the number of L1 loop iterations for the Rework matrix
            num_iters = ans_L1.nit
        
        # Do not perform L1 loop if not sufficient amount of iterating variables in equations
        else:
            
            # Set number of L1 iterations equal to 0
            num_iters = 0
            
        # Return new path values vector and number of iterations
        return self.Pvn, num_iters


# Calculate variable means and standard deviations of successful runs of all paths
class getAverages:
    
    # Initialize the class
    def __init__(self, sample, Path_vals, Run_index, mean_vals_all, std_vals_all):
        self.sample = sample
        self.Pv = Path_vals
        self.Ri = Run_index
        self.mva = mean_vals_all
        self.sva = std_vals_all
        return
    
    # Calculate variable means and standard deviations of successful runs
    def getStats(self):
        
        # Establish a counter for calculating averages and standard deviations
        count = 0
        
        # Loop through the necessary indices of each path
        for i in range(0,len(self.Ri)):
            for j in range(self.sample,len(self.Ri[i])):
                for k in range(0,len(self.Ri[i][j])):
                        
                    # Sum the values of the successful variable runs
                    self.mva[i,:] += self.Pv[i,self.Ri[i][j][k],:]
                    count += 1
                    
            # Divide the sums by the number of samples for the averages
            if count > 0:
                self.mva[i,:] = self.mva[i,:]/count
            
            # Loop back through necessary indices of each path
            for j in range(self.sample,len(self.Ri[i])):
                for k in range(0,len(self.Ri[i][j])):
                    
                    # Sum the square of the differences of the values with the averages
                    self.sva[i,:] += ((self.Pv[i,self.Ri[i][j][k],:]-self.mva[i,:])**2)
            
            # Divide the sums by the number of samples and take square root for standard deviations
            if count > 0:
                self.sva[i,:] = np.sqrt(self.sva[i,:]/count)
            
            # Reset the counter to zero
            count = 0
        
        # Return the mean values ands standard deviations for variables of each path
        return self.mva, self.sva


"""
FUNCTIONS
"""
def resultChecker(Path, bounds):
    
    # Set up empty success matrices
    Var_success = np.zeros((np.shape(Path)[0],np.shape(Path)[2]))
    Run_success = np.zeros((np.shape(Path)[0],np.shape(Path)[2]+1))
    
    # Set up empty matrix for gathering index groupings of success runs
    Run_index = [[]]*(np.shape(Path)[2]+1)
    Run_index2 = [[]]*np.shape(Path)[0]
    for i in range(0,np.shape(Path)[0]):
        Run_index2[i] = Run_index.copy()
    
    # Establish a counting variable
    count = 0
    
    # Gather success results for all of the path values
    for i in range(0,np.shape(Path)[0]):         # i loops with paths
        for j in range(0,np.shape(Path)[1]):     # j loops with runs
            for k in range(0,np.shape(Path)[2]): # k loops with variables
                if (Path[i,j,k] >= bounds[k,0]) and (Path[i,j,k] < bounds[k,1]):
                    Var_success[i,k] += 1
                    count += 1
            Run_success[i,count] += 1
            Run_index2[i][count] = Run_index2.copy()[i][count] + [j]
            count = 0
            
    # Calculate percentages from the results
    Var_success = np.around(Var_success / np.shape(Path)[1] * 100, 2)
    Run_success = np.around(Run_success / np.shape(Path)[1] * 100, 2)
    
    # Return the result percentages
    return Var_success, Run_success, Run_index2


"""
USER INPUTS
"""
# Assign number of runs for each path
runs = 10000

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

# Set tolerance for closeness of variables
tol = 5e-2

# Set initial standard deviation percentage (of the range on the bounds) for restart loop variables
std_restart = 0.1

# Set max number of analysis (L1) rework loops
l1_max = 10

# Set max number of restart (L-restart) loops
lrestart_max = 5

# Set max number of large (L4) rework loops
l4_max = 11

# Set method for the norm minimizer of the L1 rework loops
mini = 'BFGS'

# Set number of successful variables in a run to sample mean and std
sample = 9


"""
SCRIPT
"""
# Set up empty vectors and matrics
Rework_lrestart = np.zeros((np.shape(sequence)[0],runs,l4_max))
Rework_L4 = np.zeros(np.shape(sequence)[0])
mean_vals_all = np.zeros((np.shape(sequence)[0],l4_max,np.shape(bounds)[0]))
std_vals_all = np.zeros((np.shape(sequence)[0],l4_max,np.shape(bounds)[0]))
Run_success = np.zeros((np.shape(sequence)[0],np.shape(bounds)[0]+1))

# Assign input values, calculate outputs, check for conflicts, resolve
for h in range(0,l4_max):                            # h loops with L4 rework

    # Set up empty matrices
    Path_vals = np.zeros((np.shape(sequence)[0],runs,np.shape(bounds)[0]))
    Rework_L1 = np.zeros((np.shape(sequence)[0],runs,np.shape(sequence)[1]))
    
    for i in range(0,np.shape(Path_vals)[0]):        # i loops with paths
    
        # Do not loop through path if all runs are 100% successful
        if (Run_success[i,np.shape(bounds)[0]] >= 100):
            continue
    
        for j in range(0,runs):                      # j loops with runs
            
            # Define zero vectors for run
            Path_vals_new = np.zeros(np.shape(Path_vals)[2])
            mean_vals =  np.zeros(np.shape(Path_vals)[2])
            
            # Establish a sequence loop counter
            count_seq = 0
            
            # Loop over each of the analyses
            while count_seq < np.shape(sequence)[1]:
                
                # Check max restart count is not exceeded
                if Rework_lrestart[i,j,h] > lrestart_max:
                    Path_vals[i,j,:] = np.NaN
                    break
                    
                # Retrieve variables involved in analysis
                index = sequence[i][count_seq]
                variables = getVariables(analysis[index-1])
                Vars = variables.getVars()
            
                # Get random values for inputs of analysis
                random = getInput(Vars,Path_vals_new,depend[index-1][:],x)
                if (h < 1) or (np.all(mean_vals_all[i,h-1,:]) == 0):
                    if np.all(mean_vals == 0):
                        Path_vals_new = random.getUniform(bounds)
                    else:
                        Path_vals_new = random.getHybrid(bounds,mean_vals,std_restart,Rework_lrestart[i,j,h])
                else:
                    Path_vals_new = random.getNormal(mean_vals_all[i,h-1,:],std_vals_all[i,h-1,:],bounds,mean_vals,std_restart,Rework_lrestart[i,j,h])
                
                # Create function(s) for analysis with numerical inputs and variable output(s)
                func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
                expr = func.getFunc()
            
                # Evaluate analysis
                sols = sp.solve(expr)
            
                # Assign dependent variables(s) of analysis to new path values vector
                solution = assignOutput(sols,Path_vals_new,x)
                Path_vals_new = solution.solAssign()
                
                # Check for conflicts if not the first analysis in the sequence
                if count_seq > 0:
                    check = checkConflict(tol,Path_vals[i,j,:],Path_vals_new,depend[index-1][:],x)
                    conflict = check.getCheck()
            
                    # L1 loop if there are any conflicts
                    if np.any(~conflict):
                        
                        # Go straight to restart if any of the new Path values vector is NaN
                        if np.any(~np.isnan(Path_vals_new)):

                            # Gather new input values with the desired iterator
                            # Populate L1 Rework loop with the number of iterations
                            looper1 = getLoopy(Path_vals[i,j,:],Path_vals_new,Vars,analysis[index-1],depend[index-1][:],x,l1_max,mini)
                            Path_vals_new, Rework_L1[i,j,index-1] = looper1.analysisLoop()
                    
                            # Re-create function(s) for analysis with numerical inputs and variable output(s)
                            func = createFunction(analysis[index-1],Path_vals_new,Vars,depend[index-1][:],x)
                            expr = func.getFunc()
                    
                            # Re-evaluate analysis
                            sols = sp.solve(expr)
                    
                            # Re-assign dependent variable(s) of analysis to new path values vector
                            solution = assignOutput(sols,Path_vals_new,x)
                            Path_vals_new = solution.solAssign()
                    
                            # Re-check for conflicts after the L1 loops
                            check = checkConflict(tol,Path_vals[i,j,:],Path_vals_new,depend[index-1][:],x)
                            conflict = check.getCheck()
            
                        # Restart the sequence if any conflicts
                        if np.any(~conflict):
                            
                            # Gather the dependent variable value(s) for the reassignment
                            for k in range(0,np.shape(depend[index-1])[0]):
                                
                                # Assign the variable value(s) to the mean values vector
                                mean_vals[x.index(depend[index-1][k])] = Path_vals_new[x.index(depend[index-1][k])]
                        
                            # Increase the L-restart count by 1
                            Rework_lrestart[i,j,h] += 1
                            
                            # Reset the Path values vector to zeros
                            Path_vals_new = np.zeros(np.shape(Path_vals)[2])
                            
                            # Reset sequence counter to 0
                            count_seq = 0
                            
                        # Assign variable values to Path values if no conflicts
                        else:
                            # Assign a copy of new path values vector to official path values vector
                            Path_vals[i,j,:] = np.copy(Path_vals_new)
                        
                            # Increase the sequence counter by 1
                            count_seq += 1
                            
                    # Assign variable values to Path values if no conflicts
                    else:
                        # Assign a copy of new path values vector to official path values vector
                        Path_vals[i,j,:] = np.copy(Path_vals_new)
                        
                        # Increase the sequence counter by 1
                        count_seq += 1
                                
                # Do not check for conflicts if first analysis in the sequence
                else:
                    # Assign a copy of new path values vector to official path values vector
                    Path_vals[i,j,:] = np.copy(Path_vals_new)
                    
                    # Increase the sequence counter by 1
                    count_seq += 1
    
        # Increase L4 rework counter by 1
        Rework_L4[i] += 1
    
    # Calculate variable and run success
    Var_success, Run_success, Run_index = resultChecker(Path_vals,bounds)
    
    # Calculate means and standard deviations of successful runs
    averages = getAverages(sample, Path_vals, Run_index, mean_vals_all[:,h,:], std_vals_all[:,h,:])
    mean_vals_all[:,h,:], std_vals_all[:,h,:] = averages.getStats()
    
    ########################## GRAPH-RELATED CODE #############################
    
    # Graph variable success results
    fig = plt.figure(figsize=(10, 6))
    xi = np.arange(np.shape(Var_success)[1])
    x_data = ['x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9', 'x10']
    plt.title("Results for each Variable")
    plt.bar(xi-0.3, Var_success[0,:], color = 'red', width = 0.2)
    #plt.bar(xi-0.3, Var_success[0,:], color = 'darkgray', width = 0.2)
    plt.bar(xi-0.1, Var_success[1,:], color = 'green', width = 0.2)
    #plt.bar(xi-0.1, Var_success[1,:], color = 'gray', width = 0.2)
    plt.bar(xi+0.1, Var_success[2,:], color = 'blue', width = 0.2)
    #plt.bar(xi+0.1, Var_success[2,:], color = 'dimgray', width = 0.2)
    plt.bar(xi+0.3, Var_success[3,:], color = 'brown', width = 0.2)
    #plt.bar(xi+0.3, Var_success[3,:], color = 'lightgray', width = 0.2)
    plt.xticks(xi, x_data)
    plt.xlabel("Polynomial Model Variables")
    plt.ylabel("Percent of Runs within Bounds")
    plt.legend(["Path1", "Path2", "Path3", "Path4"], loc='upper right')
    plt.grid(which='major',axis='y')
    plt.show()
    
    # Graph run success results
    fig = plt.figure(figsize=(10, 6))
    xi = np.arange(np.shape(Run_success)[1])
    x_data = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    plt.title("Results for each Run")
    plt.bar(xi-0.3, Run_success[0,:], color = 'red', width = 0.2)
    #plt.bar(xi-0.3, Run_success[0,:], color = 'darkgray', width = 0.2)
    plt.bar(xi-0.1, Run_success[1,:], color = 'green', width = 0.2)
    #plt.bar(xi-0.1, Run_success[1,:], color = 'gray', width = 0.2)
    plt.bar(xi+0.1, Run_success[2,:], color = 'blue', width = 0.2)
    #plt.bar(xi+0.1, Run_success[2,:], color = 'dimgray', width = 0.2)
    plt.bar(xi+0.3, Run_success[3,:], color = 'brown', width = 0.2)
    #plt.bar(xi+0.3, Run_success[3,:], color = 'lightgray', width = 0.2)
    plt.xticks(xi, x_data)
    plt.xlabel("Number of Successful Variables in a Run")
    plt.ylabel("Percent of Runs")
    plt.legend(["Path1", "Path2", "Path3", "Path4"], loc='upper right')
    plt.grid(which='major',axis='y')
    plt.show()
    
    # Graph L1 rework results
    if h == 0:
        yi = np.zeros((np.shape(Rework_L1)[0],l4_max))
        ysums = np.zeros((np.shape(Rework_L1)[0],l4_max))
    for i in range(0,np.shape(Rework_L1)[0]):
        yscratch = np.sum(Rework_L1[i,:,:],axis=1)
        yi[i,h] = np.average(yscratch)
        if h > 0:
            ysums[i,h] = ysums[i,h-1] + yi[i,h-1]
fig, ax = plt.subplots(figsize=(10,6))
xi = ['1', '2', '3', '4']
for i in range(0,l4_max):
    if i == 0:
        ax.bar(xi,yi[:,i],width=0.4,bottom=0,label='Knowledge Loop '+str(i))
    elif np.any(yi[:,i]) != 0:
        ax.bar(xi,yi[:,i],width=0.4,bottom=ysums[:,i],label='Large Loop '+str(i))
    
ax.set_title("Analysis Loops")
ax.set_xlabel("Path")
ax.set_ylabel("Average Number of Analysis Loops")
ax.legend()
plt.grid(which='major',axis='y')
plt.show()
    
# Graph L-restart results
fig, ax = plt.subplots(figsize=(10, 6))
xi2 = ['1', '2', '3', '4']
yi2 = np.zeros(np.shape(Rework_lrestart)[0])
ysums2 = np.zeros(np.shape(Rework_lrestart)[0])
for i in range(0,np.shape(Rework_lrestart)[2]):
    for j in range(0,np.shape(Rework_lrestart)[0]):
        yi2[j] = np.average(Rework_lrestart[j,:,i])
    ax.bar(xi2,yi2,width=0.4,bottom=ysums2,label='Large Loop '+str(i))
    ysums2 = ysums2 + yi2
ax.set_title("Restart Loops")
ax.set_xlabel("Path")
ax.set_ylabel("Average Number of Restart Loops")
ax.legend()
plt.grid(which='major',axis='y')
plt.show()


# Graph L4 rework results
fig = plt.figure(figsize=(10, 6))
xi3 = ['1', '2', '3', '4']
plt.title("Knowledge Loops")
plt.bar(xi3, Rework_L4-1, color = 'red', width = 0.4)
#plt.bar(xi3, Rework_L4, color = 'darkgray', width = 0.4)
plt.yticks(np.arange(0, l4_max, step=1))
plt.xlabel("Path")
plt.ylabel("Number of Knowledge Loops")
plt.grid(which='major',axis='y')
plt.show()





















