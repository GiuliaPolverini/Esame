# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 12:06:14 2020

@author: Giulia
"""
import numpy as np
from scipy.integrate import odeint # package for ODE integration

class PendulumSimulation:
    '''
    Classe per la simulazione del sistema fisico del doppio pendolo...
    '''
    def __init__(self, state0, time_end, exp_time_sampling = 11, l = 1):
        '''
        Classe per la simulazione del sistema fisico del doppio pendolo...
        '''
        self.time = np.linspace(0.0, time_end, 2**exp_time_sampling+1)
        self.dt = self.time[1] - self.time[0]
        self.l = l
        self.state0 = state0
        self.motion = None

    def Pendulum_model(self, state, time, l):
        ''' (ESEMPIO DI DOCUMENTAZIONE)
        (descrizione) Funzione per la derivata del sistema....

        parametrs:
            state: stato iniziale del sistema
            time : campionamento temporale su cui fare la simulazione
                   ~ np.linspace(0,10, 2**9+1)
            l    : lunghezza delle sbarre del pendolo. Nell'approssimazione in cui l1 = l2.

        return:
            δa1, δa2, δp1, δp2: variazioni infinitesime delle 4 variabili del sistema.
        '''
        #TO DO: add l1, l2, m1, m2 as parametrs
        a1, a2, p1, p2 = state #a -> angle, p -> momentum
        g = 9.81

        expr1 = np.cos(a1 - a2)
        expr2 = np.sin(a1 - a2)
        expr3 = (1 + expr2**2)
        expr4 = p1 * p2 * expr2 / expr3
        expr5 = (p1**2 + 2 * p2**2 - p1 * p2 * expr1) * np.sin(2 * (a1 - a2)) / (2 * expr3**2)
        expr6 = expr4 - expr5

        δa1 = (p1 - p2 * expr1) / expr3
        δa2 = (2 * p2 - p1 * expr1) / expr3
        δp1 = (-2 * g * l * np.sin(a1) - expr6)
        δp2 = (-1 * g * l * np.sin(a2) + expr6)

        return δa1, δa2, δp1, δp2

    def simulate(self):
        self.motion = odeint(self.Pendulum_model, y0=self.state0, t=self.time, args=(self.l,))
        self.cartesian_traj()

    def cartesian_traj(self):
        '''
        Function for the conversion to cartesian coordinate of the series of values in 'self.motion'.
        '''
        a1_hat, a2_hat, p1_hat, p2_hat = self.motion.T
        self.x1 =  self.l * np.sin(a1_hat)
        self.y1 = -self.l * np.cos(a1_hat)
        self.x2 = self.x1 + self.l * np.sin(a2_hat)
        self.y2 = self.y1 - self.l * np.cos(a2_hat)