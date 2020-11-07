# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:56:19 2020

@author: Giulia
"""
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation

# class about the settings for the animation window 
class Animator:
    
    def __init__(self, pend_simulation, draw_trace = False):
        # to draw the trajectories of the masses
        self.draw_trace = draw_trace
        # lenght of the simulation (in sec) got by doing the difference between
        # two following instants in PendulumSimulation
        self.dt = pend_simulation.dt
      
        # testing the integration
        if pend_simulation.integrate is None:
            raise ValueError('Value Error: `pend_simulation` not alread integrated.')
        else:
            # saving the trajectories of the two masses
            self.x1 = pend_simulation.x1
            self.y1 = pend_simulation.y1
            self.x2 = pend_simulation.x2
            self.y2 = pend_simulation.y2 

            # boundary: fits the size of the animation window to the length of 
            # the pendulum plus its 20%
            self.bound = 0.2
            x2_extr = np.asarray([np.min(self.x2), np.max(self.x2)]) * (1 + self.bound)
            y2_extr = np.asarray([np.min(self.y2), np.max(self.y2)]) * (1 + self.bound)
           
            # set up the figure
            self.fig, self.ax = plt.subplots(figsize = (8,8))
            self.ax.set_xlim(*x2_extr) # self.ax.set_xlim(x2_extr[0], x2_extr[1])
            self.ax.set_ylim(*y2_extr)
          
            # prepare a text window for the timer
            self.time_text = self.ax.text(0.05, 0.95,
                                          'Elapsed time: {:6.2f} s'.format(0.0),
                		                  horizontalalignment = 'left',
                		                  verticalalignment = 'top',
                		                  transform = self.ax.transAxes)

            # initialize by plotting the first positions of the masses
            self.line, = self.ax.plot([0.0, self.x1[0], self.x2[0]],
                                      [0.0, self.y1[0], self.y2[0]],
                                      marker = 'o')

            # trace the whole trajectory of both the masses
            if self.draw_trace:
                self.trace1, = self.ax.plot(self.x1[0], self.y1[0])
                self.trace2, = self.ax.plot(self.x2[0], self.y2[0])
            else:
                self.trace = None

    def update(self, n):
        '''
        This function makes the text evolve: as the iteration is done for 
        regular time, the time of simulation is given by n*dt.
        
        It's a grafical function.

        Parameters
        ----------
        n : index for each number of iteration

        '''
        # fill the text window with following iterations
        # (the elapsed time must have at most six digits before comma, at most
        # two digits after comma, and it's a float number)
        self.time_text.set_text('Elapsed time: {:6.2f} s'.format(n * self.dt))

        # overwriting of each position for both the pendulums
        self.line.set_xdata([0.0, self.x1[n], self.x2[n]])
        self.line.set_ydata([0.0, self.y1[n], self.y2[n]])

        # overwriting of each trace for both the pendulums
        if self.draw_trace:
            self.trace1.set_xdata(self.x1[0: n+1])
            self.trace1.set_ydata(self.y1[0: n+1])
            self.trace2.set_xdata(self.x2[0: n+1])
            self.trace2.set_ydata(self.y2[0: n+1])
        
        return self.line, self.trace1, self.trace2
       
    def animate(self):
        '''
        This function, with FuncAnimation, generates as much frames as required
        and then passes them to the update method
        
        '''
        # cycle that iterates the sequence of frames
        self.animation = animation.FuncAnimation(fig = self.fig,
                                                 func = self.update,
                                                 frames = range(1000),
                                                 interval = 5, blit = False)