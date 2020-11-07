# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:56:19 2020

@author: Giulia
"""
import matplotlib.pylab as plt
import json

from Esame.Animator import Animator
from Esame.PendulumSimulation import PendulumSimulation

#%%----------------------------------------------------------------------------

# loading configuration  
# bring out the value I need from the config file  
with open('.\Esame\config.json') as f:
  config = json.load(f)
config['pend']['state0'] = [float(st_vl) for st_vl in config['pend']['state0']]

#creation of the simulation
pend_sim = PendulumSimulation(**config['pend']) 
pend_sim.simulate()

#creation of the animation
anim = Animator(pend_sim, **config['anim'])
anim.animate()

plt.show()

# ATTENZIONE: per clonare la mia repository togliere le parti dove c'è ESAME, 
# perché è una cosa del mio PC + attenzione agli slash --> scriverlo sul README?