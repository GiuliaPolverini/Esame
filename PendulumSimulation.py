# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 12:06:14 2020

@author: Giulia
"""
import numpy as np
from scipy.integrate import odeint # package for ODE integration

# class for the simulation of the double pendulum
class PendulumSimulation:
    
    def __init__(self, state0, time_end, exp_time_sampling = 11, l = 1):
        # testing for the length of the simulation
        if time_end <= 0:
            raise ValueError('Value Error: `time_end` must be > 0.')
        # testing for the sampling of the simulation
        elif exp_time_sampling < 0:
            raise ValueError('Value Error: `exp_time_sampling` must be >= 0.')
        else:
            self.time = np.linspace(0.0, time_end, 2**exp_time_sampling+1)

        # calculation of the length of the simulation
        self.dt = self.time[1] - self.time[0]
        self.state0 = state0
        # initialization of the motion
        self.integrate = None

        # testing of the pendulum length
        if l > 0:
            self.l = l
        else:
            raise ValueError('Value Error: l must be > 0.')

    def simulate(self):
        '''
        This function ...

        '''
        self.integrate = odeint(self._pendulum_model, y0=self.state0,
                                t=self.time, args=(self.l,))
        self.x1, self.y1, self.x2, self.y2 = self._cartesian_traj(self.integrate)

    def _pendulum_model(self, state, time, l):
        '''
        This function makes the derivatives of the system

        Parameters
        ----------
        state : initial state of the system
        time : temporal sampling su cui fare la sim
        l : length of the pendulum in the approximation where l1 = l2

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


    def _cartesian_traj(self, integrate):
        '''
        Function for the conversion to cartesian coordinate of the series of 
        values in 'self.integrate'.
        It modifies the state of the variables, building new attributes
       
        Parameters
        ----------
        integrate : 

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
    
def test_init():
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=1)
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=10000)
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=2)

def test_cartesian_traj():
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=1)

    for lenght in [200, 300, 0]:
        integrate = np.ones((lenght,4))
        a1, a2, a3, a4 = pend_sim._cartesian_traj(integrate)
        # 4variabili perché state ha 4 valori
        assert integrate.shape[0] == a1.shape[0]
