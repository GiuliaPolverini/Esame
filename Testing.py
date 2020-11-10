# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 19:50:53 2020

@author: Giulia
"""
import numpy as np
# from PendulumSimulation import PendulumSimulation
from Esame.PendulumSimulation import PendulumSimulation

def test_init():
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=1)
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=10000)
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=2)

def test_cartesian_traj():
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=1)

    for lenght in [200, 300, 0]:
        integrate = np.ones((lenght,4))
        a1, a2, a3, a4 = pend_sim._cartesian_traj(integrate)
        
        assert integrate.shape[0] == a1.shape[0]