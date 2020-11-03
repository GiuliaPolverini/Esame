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
    Qui si pone il problema
    '''
    def __init__(self, state0, time_end, exp_time_sampling = 11, l = 1):
        '''
        Classe per la simulazione del sistema fisico del doppio pendolo...
        '''
        self.time = np.linspace(0.0, time_end, 2**exp_time_sampling + 1)
        # campionamento sul tempo (=quanti istanti ci sono tra 0 e time_end;
        # faccio il campionamento chiedendo un numero che ha una potenza di 2
        # di 2 punti --> facendo la divisione per 2 ho risultati perfetti,
        # no errori macchina)
        self.dt = self.time[1] - self.time[0] # calcolo il dt per Animator
        self.l = l
        self.state0 = state0
        self.motion = None # la inizializzo vuota, la riempio sotto quando
                           # chiamo simulate
        
        # ho tutto quello che serve a Odeint per funzionare: time (time_end),
        # stato iniziale (state0) e almeno un parametro (l)

    def Pendulum_model(self, state, time, l):
        '''
        descrizione:
        Funzione per fare la derivata del sistema....
        Qui si va a risolvere il problema

        parametrs:
            state: stato iniziale del sistema
            time : campionamento temporale su cui fare la simulazione
                   ~ np.linspace(0,10, 2**9+1)
            l: lunghezza delle sbarre del pendolo. Nell'approssimazione in cui l1=l2.

        return:
            δα1, δα2, δp1, δp2: variazioni infinitesime delle 4 variabili del sistema.
        '''
        
        α1, α2, p1, p2 = state
        g = 9.81

        expr1 = np.cos(α1 - α2)
        expr2 = np.sin(α1 - α2)
        expr3 = (1 + expr2**2)
        expr4 = p1 * p2 * expr2 / expr3
        expr5 = (p1**2 + 2 * p2**2 - p1 * p2 * expr1) * np.sin(2 * (α1 - α2)) / (2 * expr3**2)
        expr6 = expr4 - expr5

        δα1 = (1 * p1 - p2 * expr1) / expr3
        δα2 = (2 * p2 - p1 * expr1) / expr3
        δp1 = (-2 * g * l * np.sin(α1) - expr6)
        δp2 = (-1 * g * l * np.sin(α2) + expr6)

        return δα1, δα2, δp1, δp2

    # metodo che risolve effettivam il moto (non avrò un moto finché non simulo)
    def simulate(self):
        self.motion = odeint(self.Pendulum_model, y0 = self.state0,
                             t = self.time, args = (self.l,))
        # risultato di OdeInt (soluzione numerica dell'integraz delle eq.diff)
        self.cartesian_traj()

    def cartesian_traj(self):
        '''
        Function for the conversion to cartesian coordinate of the series of 
        values in 'self.motion'.
        Non ha parametri in ingresso nè un return
        '''
        α1_hat, α2_hat, p1_hat, p2_hat = self.motion.T
        
        self.x1 =  self.l * np.sin(α1_hat)
        self.y1 = -self.l * np.cos(α1_hat)
        self.x2 = self.x1 + self.l * np.sin(α2_hat)
        self.y2 = self.y1 - self.l * np.cos(α2_hat)