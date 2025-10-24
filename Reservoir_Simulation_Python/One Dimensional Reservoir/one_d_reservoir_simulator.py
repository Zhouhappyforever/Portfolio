#!/usr/bin/env python
# coding: utf-8

# My Professional Code for 1D Reservoir Simulation
# 
# This is my implementation of a one-dimensional reservoir simulator using finite differences.
# It handles pressure distribution in reservoirs with variable properties and different solver methods.
# Developed by Andrew Z. Miller as part of my reservoir engineering portfolio.

# In[1]:


import numpy as np
import yaml
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt


# In[ ]:


import numpy as np
import yaml
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt

class OneDimReservoir():
    
    def __init__(self, inputs):
        '''
        My class for solving one-dimensional reservoir problems with
        finite differences. Implemented by Andrew Z. Miller.
        '''
        
        #stores input dictionary as class attribute, either read from a yaml file
        #or directly from a Python dictonary
        if isinstance(inputs, str):
            with open(inputs) as f:
                self.inputs = yaml.load(f, yaml.FullLoader)
        else:
            self.inputs = inputs
        
        #computes delta_x
        self.Nx = self.inputs['numerical']['number of grids']['x']
        self.N = self.Nx
        
        # Create arrays for properties
        self.k = self._to_array(self.inputs['reservoir']['permeability'])
        self.phi = self._to_array(self.inputs['reservoir']['porosity'])
        self.A = self._to_array(self.inputs['reservoir']['cross sectional area'])
        self.mu = self.inputs['fluid']['water']['viscosity']
        self.b_alpha = self.inputs['fluid']['water']['formation volume factor']
        self.c = self.inputs['fluid']['water']['compressibility']
        
        if 'delta x' in self.inputs['numerical']:
            self.dx = self._to_array(self.inputs['numerical']['delta x'])
        else:
            self.dx = np.full(self.N, self.inputs['reservoir']['length'] / self.N)
        
        #gets delta_t from inputs
        self.delta_t = self.inputs['numerical']['time step']
        
        #applies the initial reservoir pressues to self.p
        self.apply_initial_conditions()
        
        #calls fill matrix method (must be completely implemented to work)
        self.fill_matrices()
        
        #create an empty list for storing data if plots are requested
        if 'plots' in self.inputs:
            self.p_plot = []
        
        return
    
    def _to_array(self, val):
        if isinstance(val, (list, tuple, np.ndarray)):
            arr = np.array(val)
            if len(arr) != self.N:
                raise ValueError("Array length must match number of grids")
            return arr
        else:
            return np.full(self.N, val)
    
    def compute_transmissibility(self,i,j):
        '''
        Computes the transmissibility. Part of my reservoir simulation code.
        '''
        if i == j:
            return self.k[i] * self.A[i] / (self.mu * self.b_alpha * self.dx[i])
        else:
            trans = 2 * self.k[i] * self.k[j] * self.A[i] * self.A[j] / (self.k[i] * self.A[i] * self.dx[j] + self.k[j] * self.A[j] * self.dx[i])
            return trans / (self.mu * self.b_alpha)
    
    
    def compute_accumulation(self, i):
        '''
        Computes the accumulation. Part of my reservoir simulation code.
        '''
        Vb = self.A[i] * self.dx[i]
        acc = self.phi[i] * self.c * Vb / self.b_alpha
        return acc
    
    def fill_matrices(self):
        '''
        Fills the matrices T, B, and \vec{Q} and applies boundary
        conditions. Part of my reservoir simulation code.
        '''
        conv = self.inputs['conversion factor']
        
        N = self.N
        
        acc = [self.compute_accumulation(i) for i in range(N)]
        self.B = scipy.sparse.diags([acc], [0], format='csr')
        
        self.Q = np.zeros(N)
        
        # Compute interblock transmissibilities
        if N > 1:
            T_half = np.array([self.compute_transmissibility(i, i+1) * conv for i in range(N-1)])
        else:
            T_half = np.array([])
        
        # Build diagonals for T
        if N > 1:
            lower = -T_half
            upper = -T_half
            main_diag = np.zeros(N)
            main_diag[0] = T_half[0]
            main_diag[-1] = T_half[-1]
            for ii in range(1, N-1):
                main_diag[ii] = T_half[ii-1] + T_half[ii]
        else:
            lower = np.array([])
            upper = np.array([])
            main_diag = np.zeros(N)
        
        self.T = scipy.sparse.diags([lower, main_diag, upper], [-1, 0, 1], format='csr')
        
        # Convert to lil for modifications
        self.T = self.T.tolil()
        
        # Apply boundary conditions
        left_bc = self.inputs['boundary conditions']['left']
        if left_bc['type'] == 'prescribed pressure':
            p_left = left_bc['value']
            T_left = 2 * self.compute_transmissibility(0, 0) * conv
            self.T[0, 0] += T_left
            self.Q[0] += T_left * p_left
        elif left_bc['type'] == 'prescribed flux':
            q_left = left_bc['value']
            self.Q[0] += q_left
        
        right_bc = self.inputs['boundary conditions']['right']
        if right_bc['type'] == 'prescribed pressure':
            p_right = right_bc['value']
            T_right = 2 * self.compute_transmissibility(self.N-1, self.N-1) * conv
            self.T[-1, -1] += T_right
            self.Q[-1] += T_right * p_right
        elif right_bc['type'] == 'prescribed flux':
            q_right = right_bc['value']
            self.Q[-1] += q_right
        
        # Convert back to csr
        self.T = self.T.tocsr()
        
        return
    
    
    def apply_initial_conditions(self):
        '''
        Applies initial pressures to self.p. Part of my reservoir simulation code.
        '''
        
        N = self.N
        
        self.p = np.full(N, self.inputs['initial conditions']['pressure'])
        
        return
    
    
    def solve_one_step(self):
        '''
        Solve one time step using either the implicit or explicit method.
        Part of my reservoir simulation code.
        '''
        solver_input = self.inputs['numerical']['solver']
        if isinstance(solver_input, str):
            method = solver_input
            theta = None
        elif isinstance(solver_input, dict):
            method = 'mixed method'
            theta = solver_input['mixed method']['theta']
        else:
            raise ValueError("Invalid solver input")
        
        dt = self.delta_t
        p = self.p
        B = self.B
        T = self.T
        Q = self.Q

        if method == 'explicit':
            self.p = p + dt * np.linalg.solve(B.toarray(), Q - T.dot(p))
        elif method == 'implicit':
            self.p = np.linalg.solve((B.toarray() / dt) + T.toarray(), (B.toarray() / dt).dot(p) + Q)
        elif method == 'mixed method':
            self.p = np.linalg.solve((B.toarray() / dt) + (1 - theta) * T.toarray(), (B.toarray() / dt - theta * T.toarray()).dot(p) + Q)

        
        return
    
    
    def solve(self):
        '''
        Solves until "number of time steps". Part of my reservoir simulation code.
        '''
        
        for i in range(self.inputs['numerical']['number of time steps']):
            self.solve_one_step()
            
            if 'plots' in self.inputs and i % self.inputs['plots']['frequency'] == 0:
                self.p_plot += [self.get_solution().copy()]
        
        return
    
    def plot(self):
        '''
        Crude plotting function. Plots pressure as a function of grid block #.
        Part of my reservoir simulation code.
        '''
        
        if self.p_plot is not None:
            for i in range(len(self.p_plot)):
                plt.plot(self.p_plot[i])
        
        return
    
    def get_solution(self):
        '''
        Returns solution vector. Part of my reservoir simulation code.
        '''
        return self.p


# # Example code execution
# 
# If you'd like to run the code in the notebook, perhaps creating a crude plot of the output, you can uncomment the following lines of code in the cell below.  You can also inspect the contents of `inputs.yml` and change the parameters to see how the solution is affected.

# In[4]:


#import matplotlib.pyplot as plt
#implicit = OneDimReservoir('inputs.yml')
#implicit.solve()
#implicit.plot()