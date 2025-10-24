#!/usr/bin/env python
# coding: utf-8

# Here is my code for reservoir simulation in a professional version.
# This script implements the Buckley-Leverett semi-analytical solution for two-phase flow.
# It computes relative permeabilities, fractional flow, and saturation profiles with shock front correction.

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import yaml

class BuckleyLeverett:
    """
    Professional implementation of the Buckley-Leverett semi-analytical solution.
    
    This class models water-oil displacement in a 1D reservoir, computing relative permeabilities
    using the Corey-Brooks model, fractional flow functions, and the shock front saturation via
    nonlinear optimization. It provides methods for plotting fractional flow and saturation profiles.
    
    Attributes:
        inputs (dict): Input parameters for the model.
        Sor (float): Residual oil saturation.
        Swc (float): Critical water saturation.
        nw (float): Water Corey-Brooks exponent.
        no (float): Oil Corey-Brooks exponent.
        krw_max (float): Maximum water relative permeability.
        kro_max (float): Maximum oil relative permeability.
        mu_o (float): Oil viscosity.
        mu_w (float): Water viscosity.
        Swi (float): Initial water saturation.
        step (float): Saturation step size for computations.
    """
    
    def __init__(self, inputs):
        """
        Initializes the model with input parameters.
        
        Args:
            inputs (str or dict): Path to YAML file or dictionary of input parameters.
        """
        if isinstance(inputs, str):
            with open(inputs) as f:
                self.inputs = yaml.safe_load(f)
        else:
            self.inputs = inputs
            
        self.Sor = self.inputs['reservoir']['oil']['residual saturation']
        self.Swc = self.inputs['reservoir']['water']['critical saturation']
        self.nw = self.inputs['reservoir']['water']['corey-brooks exponent']
        self.no = self.inputs['reservoir']['oil']['corey-brooks exponent']
        self.krw_max = self.inputs['reservoir']['water']['max relative permeability']
        self.kro_max = self.inputs['reservoir']['oil']['max relative permeability']
        
        self.mu_o = self.inputs['fluid']['oil']['viscosity']
        self.mu_w = self.inputs['fluid']['water']['viscosity']
        
        self.Swi = self.inputs['initial conditions']['water saturation']
        
        self.step = 0.01
    
    def water_rel_perm(self, S):
        """Computes water relative permeability using Corey-Brooks model."""
        return self.krw_max * ((S - self.Swc) / (1 - self.Sor - self.Swc)) ** self.nw
    
    def oil_rel_perm(self, S):
        """Computes oil relative permeability using Corey-Brooks model."""
        return self.kro_max * ((1 - S - self.Sor) / (1 - self.Sor - self.Swc)) ** self.no
    
    def fractional_flow(self, S):
        """Computes fractional flow of water."""
        krw = self.water_rel_perm(S)
        kro = self.oil_rel_perm(S)
        return (krw / self.mu_w) / (krw / self.mu_w + kro / self.mu_o)

    def d_fractional_flow_dS(self, S):
        """Computes derivative of fractional flow with respect to saturation."""
        a = 1 - self.Sor - self.Swc
        Se_w = (S - self.Swc) / a
        Se_o = (1 - S - self.Sor) / a
        
        krw = self.krw_max * Se_w ** self.nw
        kro = self.kro_max * Se_o ** self.no
        
        lw = krw / self.mu_w
        lo = kro / self.mu_o
        
        dkrw_dS = self.krw_max * self.nw * Se_w ** (self.nw - 1) * (1 / a)
        dkro_dS = self.kro_max * self.no * Se_o ** (self.no - 1) * (-1 / a)
        
        dlw_dS = dkrw_dS / self.mu_w
        dlo_dS = dkro_dS / self.mu_o
        
        return (dlw_dS * lo - lw * dlo_dS) / (lw + lo)**2
    
    def compute_saturation_front(self):
        """Solves for saturation at the shock front using nonlinear optimization."""
        def g(S):
            f_wf = self.fractional_flow(S)
            return f_wf / (S - self.Swi) - self.d_fractional_flow_dS(S)
        
        x_mid = (self.Swi + (1 - self.Sor)) / 2
        return scipy.optimize.root_scalar(g, method='secant', x0=self.Swi + 0.1, x1=x_mid).root
        
    def compute_saturation_profile(self):
        """Computes the full (inadmissible) saturation profile."""
        S = np.arange(self.Swi + self.step, (1 - self.Swc), self.step)
        x = self.d_fractional_flow_dS(S)
        return (x, S)
    
    def plot_fractional_flow(self):
        """Plots the fractional flow curve."""
        S = np.arange(self.Swi + self.step, (1-self.Swc), self.step)
        f = self.fractional_flow(S)
        plt.plot(S, f)
        plt.xlabel('$S_w$')
        plt.ylabel('$f$')
        plt.title('Fractional Flow Curve')
        plt.show()
    
    def plot_full_saturation_profile(self):
        """Plots the full saturation profile without shock correction."""
        x, S = self.compute_saturation_profile()
        plt.plot(x, S)
        plt.ylabel('$S_w$')
        plt.xlabel('$x$')
        plt.title('Full Saturation Profile')
        plt.show()
        
    def plot_saturation_profile(self, t):
        """Plots the corrected saturation profile at time t with shock front."""
        x, S = self.compute_saturation_profile()
        Swf = self.compute_saturation_front()
        S1 = S[S > Swf]
        x1 = x[S > Swf] * t
        xD = self.d_fractional_flow_dS(Swf) * t
        S_plot = np.concatenate((S1[::-1], np.array([Swf, self.Swi]), np.array([self.Swi])))
        x_plot = np.concatenate((x1[::-1], np.array([xD, xD]), np.array([1.0])))
        plt.plot(x_plot, S_plot)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$S_w$')
        plt.title(f'Saturation Profile at t={t}')
        plt.show()