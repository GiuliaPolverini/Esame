# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 12:06:14 2020

@author: Giulia
"""
import numpy as np
from scipy.integrate import odeint # package for ODE integration

# class for the simulation of the double pendulum
class PendulumSimulation:
    
    '''
    Simulation of the physical system of a double pendulum throught the integration
    of a set of four first order differential equations, two for each degree of freedom.
    
    '''
    
    def __init__(self, state0, time_end, exp_time_sampling = 11, l = 1):
        
        '''
        Parameters
        ----------
        state0 : initial state of the system
        time_end : 
        exp_time_sampling : 
        l : 

        Returns
        -------
        None.

        '''
        
        if time_end <= 0:
            raise ValueError('Value Error: `time_end` must be > 0.')
        elif exp_time_sampling < 0:
            raise ValueError('Value Error: `exp_time_sampling` must be >= 0.')
        else:
            self.time = np.linspace(0.0, time_end, 2**exp_time_sampling+1)

        # calculation of the length of the simulation
        self.dt = self.time[1] - self.time[0]
        self.state0 = state0
       
        self.integrate = None

        if l > 0:
            self.l = l
        else:
            raise ValueError('Value Error: l must be > 0.')

    def simulate(self):
        
        '''
        Numerical integration of the motion of the pendulum according to the
        physical model, starting by certain initial condition 'state0', over 
        a certain interval 'time'.

        '''
        
        self.integrate = odeint(self.pendulum_model, y0=self.state0,
                                t=self.time, args=(self.l,))
        self.x1, self.y1, self.x2, self.y2 = self.cartesian_traj(self.integrate)

    def pendulum_model(self, state, time, l):
        
        '''
        This function makes the derivatives of the system.

        Parameters
        ----------
        state : initial state of the system
        time : temporal sampling of the simulation
        l : length of the pendulum, in the approximation where l1 = l2

        Returns
        -------
        δθ1, δθ2, δp1, δp2 : infinitesimal variations of the 4 variables of the system
        
        '''

        θ1, θ2, p1, p2 = state 
        g = 9.81

        expr1 = np.cos(θ1 - θ2)
        expr2 = np.sin(θ1 - θ2)
        expr3 = (1 + expr2**2)
        expr4 = p1 * p2 * expr2 / expr3
        expr5 = (p1**2 + 2 * p2**2 - p1 * p2 * expr1) * np.sin(2 * (θ1 - θ2)) / (2 * expr3**2)
        expr6 = expr4 - expr5

        δθ1 = (1 * p1 - p2 * expr1) / expr3
        δθ2 = (2 * p2 - p1 * expr1) / expr3
        δp1 = (-2 * g * l * np.sin(θ1) - expr6)
        δp2 = (-1 * g * l * np.sin(θ2) + expr6)

        return δθ1, δθ2, δp1, δp2


    def cartesian_traj(self, integrate):
        '''
        Function for the conversion from angular to cartesian coordinates of 
        the series of values in 'self.integrate'.
       
        Parameters
        ----------
        integrate : temporal series of two angular variables and their momentums

        Returns
        -------
        x1, y1, x2, y2 : converted variables
        
        '''
        
        θ1_hat, θ2_hat, p1_hat, p2_hat = integrate.T
        x1 =  self.l * np.sin(θ1_hat)
        y1 = -self.l * np.cos(θ1_hat)
        x2 = x1 + self.l * np.sin(θ2_hat)
        y2 = y1 - self.l * np.cos(θ2_hat)

        return x1, y1, x2, y2 
