# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:56:19 2020

@author: Giulia
"""
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation

class Animator:
    
    def __init__(self, pend_simulation, draw_trace = False):
        self.draw_trace = draw_trace
        self.dt = pend_simulation.dt

        if pend_simulation.motion is None:
            print('Moto del pendolo non simulato.')
        else:
            # Saving the trajectories
            self.x1 = pend_simulation.x1
            self.y1 = pend_simulation.y1
            self.x2 = pend_simulation.x2
            self.y2 = pend_simulation.y2

            self.bound = 0.2 #boundary percentage

            x2_extr = np.asarray([np.min(self.x2), np.max(self.x2)]) * (1 + self.bound)
            y2_extr = np.asarray([np.min(self.y2), np.max(self.y2)]) * (1 + self.bound)

            # set up the figure
            self.fig, self.ax = plt.subplots( figsize=(8,8))
            self.ax.set_xlim(*x2_extr) #<- self.ax.set_xlim(x2_extr[0], x2_extr[1])
            self.ax.set_ylim(*y2_extr)

            # prepare a text window for the timer
            self.time_text = self.ax.text(0.05, 0.95, 'Elapsed time: {:6.2f} s'.format(0.0),
                		     horizontalalignment = 'left',
                		     verticalalignment = 'top',
                		     transform = self.ax.transAxes)

            # initialize by plotting the last position of the trajectory
            self.line, = self.ax.plot([0.0, self.x1[0], self.x2[0]],
                                      [0.0, self.y1[0], self.y2[0]],
                                      marker = 'o')

            # trace the whole trajectory of the second pendulum mass
            if self.draw_trace:
                self.trace, = self.ax.plot(self.x2[0], self.y2[0])
            else:
                self.trace = None

    def update(self, n):
        self.time_text.set_text('Elapsed time: {:6.2f} s'.format( n * self.dt))

        # self.time_text.set_text('Frame: {} '.format(n))
        self.line.set_xdata([0.0, self.x1[n], self.x2[n]])
        self.line.set_ydata([0.0, self.y1[n], self.y2[n]])

        if self.draw_trace:
            self.trace.set_xdata(self.x2[0: n+1])
            self.trace.set_ydata(self.y2[0: n+1])
        return self.line, self.trace

    def animate(self):
        self.animation = animation.FuncAnimation(fig    = self.fig,
                                                 func   = self.update,
                                                 frames = range(1000),
                                                 interval = 5, blit = False)