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
# estraggo i valori che mi servono dal file di configurazione e li passo sotto
with open('.\Esame\config.json') as f:
  config = json.load(f)
config['pend']['state0'] = [float(st_vl) for st_vl in config['pend']['state0']]

#pendulum simulation <-- creazione della simulazione
pend_sim = PendulumSimulation(**config['pend']) # "key = value" per tutti gli 
                                                # elementi del dizionario
pend_sim.simulate() # simulazione vera e propria; simulate è una funzione che
                    # cambia lo stato dell'oggetto pend_sim, non ritorna niente

#pendulum animation
anim = Animator(pend_sim, **config['anim'])
anim.animate()
plt.show()

# ATTENZIONE: per clonare la mia repository togliere le parti dove c'è ESAME, perché è una
# cosa del mio PC + attenzione agli slash --> scriverlo sul README