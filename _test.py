# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 19:50:53 2020

@author: Giulia
"""
import numpy as np
from PendulumSimulation import PendulumSimulation
#from Esame.PendulumSimulation import PendulumSimulation

def test_init():
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=1)
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=1000)
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=2)

def test_down_pendulum():
    #test for evolution of the pendulum in its equilibrium position
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=10000)
    pend_sim.simulate()
    assert np.unique(pend_sim.integrate) == np.array([0])
    l = pend_sim.l
    for cor, val in zip([pend_sim.x1, pend_sim.y1, pend_sim.x2, pend_sim.y2],(0,-l, 0, -2*l)):
        assert np.unique(cor) == np.array([val])

def test_up_pendulum():
    #test for evolution of the pendulum in its equilibrium position
    pend_sim = PendulumSimulation(state0 = (np.pi,np.pi,0,0), time_end=20, exp_time_sampling=13)
    pend_sim.simulate()
    assert np.unique(np.around(pend_sim.integrate[0:1000,2:4], 5 )) == np.array([0])

def test_fold_pendulum_1():
    #test for evolution of the pendulum in its equilibrium position
    pend_sim = PendulumSimulation(state0 = (0,np.pi,0,0), time_end=20, exp_time_sampling=13)
    pend_sim.simulate()
    assert np.unique(np.around(pend_sim.integrate[0:1000,2:4], 5 )) == np.array([0])

def test_fold_pendulum_2():
    pend_sim = PendulumSimulation(state0 = (np.pi,0,0,0), time_end=20, exp_time_sampling=13)
    pend_sim.simulate()
    assert np.unique(np.around(pend_sim.integrate[0:1000,2:4], 5 )) == np.array([0])

def test_simulate():
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=10000)
    assert pend_sim.integrate == None
    pend_sim.simulate()
    assert pend_sim.integrate is not None

    for var in [pend_sim.x1, pend_sim.y1, pend_sim.x2, pend_sim.y2]:
        assert var is not None

def test_cartesian_traj():
    pend_sim = PendulumSimulation(state0 = (0,0,0,0), time_end=1)

    for lenght in np.random.randint(low=0, high=10_000, size=50):
        integrate = np.ones((lenght,4))
        assert integrate.shape ==  np.asarray(pend_sim.cartesian_traj(integrate)).T.shape
