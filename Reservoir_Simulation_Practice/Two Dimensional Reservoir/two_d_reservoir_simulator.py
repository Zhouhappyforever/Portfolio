#!/usr/bin/env python
# coding: utf-8

# My Professional Code for 2D Reservoir Simulation
# 
# This is my implementation of a two-dimensional reservoir simulator using finite differences.
# It handles pressure distribution in 2D reservoirs with variable properties, boundary conditions, and wells.
# Developed by Andrew Z. Miller as part of my reservoir engineering portfolio.

import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import yaml

class TwoDimReservoir():
    
    def __init__(self, inputs):
        '''
        My class for solving two-dimensional reservoir problems with
        finite differences. Implemented by Andrew Z. Miller.
        '''
        if isinstance(inputs, str):
            with open(inputs) as f:
                self.inputs = yaml.load(f, yaml.FullLoader)
        else:
            self.inputs = inputs
        
        self.Nx = self.inputs['numerical']['number of grids']['x']
        self.Ny = self.inputs['numerical']['number of grids']['y']
        self.N = self.Nx * self.Ny
        self.L = self.inputs['reservoir']['length']
        self.H = self.inputs['reservoir']['height']
        
        # Delta_x and delta_y as arrays
        if 'delta x' in self.inputs['numerical']:
            val = self.inputs['numerical']['delta x']
            self.delta_x = np.array(val) if isinstance(val, list) else np.full(self.Nx, val)
        else:
            self.delta_x = np.full(self.Nx, self.L / self.Nx)
        
        if 'delta y' in self.inputs['numerical']:
            val = self.inputs['numerical']['delta y']
            self.delta_y = np.array(val) if isinstance(val, list) else np.full(self.Ny, val)
        else:
            self.delta_y = np.full(self.Ny, self.H / self.Ny)
        
        # Properties as arrays
        self.k = self._get_property('reservoir', 'permeability')
        self.phi = self._get_property('reservoir', 'porosity')
        self.depth = self.inputs['reservoir']['depth']
        
        self.delta_t = self.inputs['numerical']['time step']
        self.mu = self.inputs['fluid']['water']['viscosity']
        self.b_alpha = self.inputs['fluid']['water']['formation volume factor']
        self.c = self.inputs['fluid']['water']['compressibility']
        
        self.apply_initial_conditions()
        
        self.fill_matrices()
        
        if 'plots' in self.inputs:
            self.p_plot = []
        
        return
    
    def _get_property(self, section, key):
        val = self.inputs[section][key]
        if isinstance(val, str):
            # Read from file if string
            with open(val, 'r') as f:
                data = [float(line.strip()) for line in f if line.strip()]
            if len(data) != self.N:
                raise ValueError(f"Data in {val} does not match number of grids")
            return np.array(data)
        elif isinstance(val, (list, tuple, np.ndarray)):
            arr = np.array(val)
            if len(arr) != self.N:
                raise ValueError(f"Array length for {key} must match number of grids")
            return arr
        else:
            return np.full(self.N, val)
    
    def get_block_number(self, x, y):
        ix = np.searchsorted(np.cumsum(self.delta_x), x, side='right') - 1
        iy = np.searchsorted(np.cumsum(self.delta_y), y, side='right') - 1
        if ix < 0 or ix >= self.Nx or iy < 0 or iy >= self.Ny:
            raise ValueError("Location out of bounds")
        return iy * self.Nx + ix
    
    def get_coords(self, block):
        iy = block // self.Nx
        ix = block % self.Nx
        return ix, iy
    
    def compute_transmissibility(self, i, j):
        '''
        Computes the transmissibility between blocks i and j. Part of my reservoir simulation code.
        '''
        if i == j:
            return 0  # Self-transmissibility not used
        # Determine if horizontal or vertical
        ix_i, iy_i = self.get_coords(i)
        ix_j, iy_j = self.get_coords(j)
        if iy_i == iy_j:  # Horizontal
            k_i = self.k[i]
            k_j = self.k[j]
            dy = self.delta_y[iy_i]
            dx_i = self.delta_x[ix_i]
            dx_j = self.delta_x[ix_j]
            A = dy * self.depth
            trans = 2 * k_i * k_j * A / (k_i * dx_j + k_j * dx_i)
        elif ix_i == ix_j:  # Vertical
            k_i = self.k[i]
            k_j = self.k[j]
            dx = self.delta_x[ix_i]
            dy_i = self.delta_y[iy_i]
            dy_j = self.delta_y[iy_j]
            A = dx * self.depth
            trans = 2 * k_i * k_j * A / (k_i * dy_j + k_j * dy_i)
        else:
            return 0
        return trans / (self.mu * self.b_alpha)
    
    def compute_accumulation(self, i):
        '''
        Computes the accumulation for block i. Part of my reservoir simulation code.
        '''
        ix, iy = self.get_coords(i)
        Vb = self.delta_x[ix] * self.delta_y[iy] * self.depth
        acc = self.phi[i] * self.c * Vb / self.b_alpha
        return acc
    
    def fill_matrices(self):
        '''
        Fills the matrices T, B, and Q and applies boundary conditions and wells. Part of my reservoir simulation code.
        '''
        conv = self.inputs['conversion factor']
        
        acc = [self.compute_accumulation(i) for i in range(self.N)]
        self.B = scipy.sparse.diags([acc], [0], format='csr')
        
        self.Q = np.zeros(self.N)
        
        # Build T matrix
        row = []
        col = []
        data = []
        
        for i in range(self.N):
            ix, iy = self.get_coords(i)
            
            # Left neighbor
            if ix > 0:
                j = i - 1
                T_ij = self.compute_transmissibility(i, j) * conv
                row.extend([i, i])
                col.extend([j, i])
                data.extend([-T_ij, T_ij])
            # Right neighbor
            if ix < self.Nx - 1:
                j = i + 1
                T_ij = self.compute_transmissibility(i, j) * conv
                row.extend([i, i])
                col.extend([j, i])
                data.extend([-T_ij, T_ij])
            # Bottom neighbor
            if iy > 0:
                j = i - self.Nx
                T_ij = self.compute_transmissibility(i, j) * conv
                row.extend([i, i])
                col.extend([j, i])
                data.extend([-T_ij, T_ij])
            # Top neighbor
            if iy < self.Ny - 1:
                j = i + self.Nx
                T_ij = self.compute_transmissibility(i, j) * conv
                row.extend([i, i])
                col.extend([j, i])
                data.extend([-T_ij, T_ij])
        
        self.T = scipy.sparse.coo_matrix((data, (row, col)), shape=(self.N, self.N)).tolil()
        
        # Boundary conditions
        total_area_x = self.depth * self.H
        total_area_y = self.depth * self.L
        
        left_bc = self.inputs['boundary conditions']['left']
        if left_bc['type'] == 'prescribed pressure':
            p_left = left_bc['value']
            for iy in range(self.Ny):
                i = iy * self.Nx
                dy = self.delta_y[iy]
                area = dy * self.depth
                T_left = 2 * self.k[i] * area / (self.mu * self.b_alpha * self.delta_x[0]) * conv
                self.T[i, i] += T_left
                self.Q[i] += T_left * p_left
        elif left_bc['type'] = 'prescribed flux':
            for iy in range(self.Ny):
                i = iy * self.Nx
                dy = self.delta_y[iy]
                area = dy * self.depth
                q_block = left_bc['value'] * area / total_area_x
                self.Q[i] += q_block
        
        right_bc = self.inputs['boundary conditions']['right']
        if right_bc['type'] == 'prescribed pressure':
            p_right = right_bc['value']
            for iy in range(self.Ny):
                i = iy * self.Nx + self.Nx - 1
                dy = self.delta_y[iy]
                area = dy * self.depth
                T_right = 2 * self.k[i] * area / (self.mu * self.b_alpha * self.delta_x[-1]) * conv
                self.T[i, i] += T_right
                self.Q[i] += T_right * p_right
        elif right_bc['type'] == 'prescribed flux':
            for iy in range(self.Ny):
                i = iy * self.Nx + self.Nx - 1
                dy = self.delta_y[iy]
                area = dy * self.depth
                q_block = right_bc['value'] * area / total_area_x
                self.Q[i] += q_block
        
        bottom_bc = self.inputs['boundary conditions']['bottom']
        if bottom_bc['type'] == 'prescribed pressure':
            p_bottom = bottom_bc['value']
            for ix in range(self.Nx):
                i = ix
                dx = self.delta_x[ix]
                area = dx * self.depth
                T_bottom = 2 * self.k[i] * area / (self.mu * self.b_alpha * self.delta_y[0]) * conv
                self.T[i, i] += T_bottom
                self.Q[i] += T_bottom * p_bottom
        elif bottom_bc['type'] == 'prescribed flux':
            for ix in range(self.Nx):
                i = ix
                dx = self.delta_x[ix]
                area = dx * self.depth
                q_block = bottom_bc['value'] * area / total_area_y
                self.Q[i] += q_block
        
        top_bc = self.inputs['boundary conditions']['top']
        if top_bc['type'] == 'prescribed pressure':
            p_top = top_bc['value']
            for ix in range(self.Nx):
                i = ix + (self.Ny - 1) * self.Nx
                dx = self.delta_x[ix]
                area = dx * self.depth
                T_top = 2 * self.k[i] * area / (self.mu * self.b_alpha * self.delta_y[-1]) * conv
                self.T[i, i] += T_top
                self.Q[i] += T_top * p_top
        elif top_bc['type'] == 'prescribed flux':
            for ix in range(self.Nx):
                i = ix + (self.Ny - 1) * self.Nx
                dx = self.delta_x[ix]
                area = dx * self.depth
                q_block = top_bc['value'] * area / total_area_y
                self.Q[i] += q_block

        # Wells
        if 'wells' in self.inputs:
            wells = self.inputs['wells']
            if 'rate' in wells:
                rate_locs = wells['rate']['locations']
                rate_vals = wells['rate']['values']
                for loc, q in zip(rate_locs, rate_vals):
                    block = self.get_block_number(loc[0], loc[1])
                    self.Q[block] += q
            if 'bhp' in wells:
                bhp_locs = wells['bhp']['locations']
                bhp_vals = wells['bhp']['values']
                bhp_radii = wells['bhp']['radii']
                skin = wells['bhp'].get('skin factor', 0.0)
                for loc, pw, rw in zip(bhp_locs, bhp_vals, bhp_radii):
                    block = self.get_block_number(loc[0], loc[1])
                    ix, iy = self.get_coords(block)
                    dx = self.delta_x[ix]
                    dy = self.delta_y[iy]
                    kx = self.k[block]
                    ky = kx  # assuming isotropic
                    req = 0.28 * np.sqrt(np.sqrt(ky / kx) * dx**2 + np.sqrt(kx / ky) * dy**2) / ((ky / kx)**0.25 + (kx / ky)**0.25)
                    ln_req_rw = np.log(req / rw) + skin
                    J = conv * 2 * np.pi * np.sqrt(kx * ky) * self.depth / (self.mu * ln_req_rw)
                    self.T[block, block] += J
                    self.Q[block] += J * pw
        
        self.T = self.T.tocsr()
        return
    
    def apply_initial_conditions(self):
        '''
        Applies initial pressures to self.p. Part of my reservoir simulation code.
        '''
        self.p = np.full(self.N, self.inputs['initial conditions']['pressure'])
        return
    
    def solve_one_step(self):
        '''
        Solve one time step using the specified method. Part of my reservoir simulation code.
        '''
        solver_input = self.inputs['numerical']['solver']
        dt = self.delta_t
        p_old = self.p.copy()
        
        if isinstance(solver_input, str):
            if solver_input == 'explicit':
                theta = 1.0
            elif solver_input == 'implicit':
                theta = 0.0
            else:
                raise ValueError(f"Unknown solver: {solver_input}")
        elif isinstance(solver_input, dict) and 'mixed method' in solver_input:
            theta = solver_input['mixed method']['theta']
        else:
            raise ValueError("Invalid solver input")
        
        lhs = self.B / dt + (1 - theta) * self.T
        rhs = (self.B / dt) @ p_old + theta * self.T @ p_old + self.Q
        self.p = scipy.sparse.linalg.spsolve(lhs, rhs)
        
        return
    
    def solve(self):
        '''
        Solves for the specified number of time steps. Part of my reservoir simulation code.
        '''
        n_steps = self.inputs['numerical']['number of time steps']
        frequency = self.inputs['plots']['frequency'] if 'plots' in self.inputs else n_steps + 1
        self.p_plot = []
        for step in range(n_steps):
            self.solve_one_step()
            if (step + 1) % frequency == 0:
                self.p_plot.append(self.p.copy())
        return
    
    def get_solution(self):
        '''
        Returns the current solution vector. Part of my reservoir simulation code.
        '''
        return self.p
    
    def plot(self):
        '''
        Plots the pressure contour. Part of my reservoir simulation code.
        '''
        if not self.p_plot:
            p_2d = self.p.reshape(self.Ny, self.Nx)
        else:
            p_2d = self.p_plot[-1].reshape(self.Ny, self.Nx)
        k_2d = self.k.reshape(self.Ny, self.Nx)
        p_masked = np.ma.masked_where(k_2d == 0, p_2d)
        
        x_edges = np.cumsum(np.insert(self.delta_x, 0, 0))
        y_edges = np.cumsum(np.insert(self.delta_y, 0, 0))
        x_centers = (x_edges[:-1] + x_edges[1:]) / 2
        y_centers = (y_edges[:-1] + y_edges[1:]) / 2
        
        X, Y = np.meshgrid(x_centers, y_centers)
        
        plt.figure(figsize=(8, 6))
        cont = plt.contourf(X, Y, p_masked, cmap='viridis', levels=50)
        plt.colorbar(cont)
        plt.xlabel('x (ft)')
        plt.ylabel('y (ft)')
        plt.title('Pressure Contour in Reservoir')
        plt.show()