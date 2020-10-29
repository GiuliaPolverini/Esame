# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:56:19 2020

@author: Giulia
"""
import matplotlib.pylab as plt
import json

from Esame.Animator import Animator
from Esame.PendulumSimulation import PendulumSimulation

#%%-----------------------------------------------------------------------------

#loading configuration    
with open('.\Esame\config.json') as f:
  config = json.load(f)
config['pend']['state0'] = [float(st_vl) for st_vl in config['pend']['state0']]

#pendulum simulation
pend_sim = PendulumSimulation(**config['pend']) #<- "key = value" per tutti gli elementi del dizionario
pend_sim.simulate()

#pendulum animation
anim = Animator(pend_sim, **config['anim'])
anim.animate()
plt.show()
