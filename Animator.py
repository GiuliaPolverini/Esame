# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:56:19 2020

@author: Giulia
"""
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation


# class to create the animation
class Animator:
    
    def __init__(self, pend_simulation, n_frames, draw_trace = False):
        
        '''
        Parameters
        ----------
        pend_simulation : simulation to animate
        n_frames : integer variable which represents the steps of the animation
        draw_trace : boolean variable to draw the trajectory of the masses

        Returns
        -------
        None.

        '''
        self.n_frames = n_frames
        self.draw_trace = draw_trace
        # size of the temporal sampling
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
            self.ax.set_xlim(*x2_extr)
            self.ax.set_ylim(*y2_extr)
          
            # prepare a text window for the timer
            self.time_text = self.ax.text(0.05, 0.95,'Elapsed time: {:6.2f} s'.format(0.0),
                		                  horizontalalignment = 'left',
                		                  verticalalignment = 'top',
                		                  transform = self.ax.transAxes)

            # the body of the pendulum is initialized by plotting the first 
            # positions of the masses
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
        This function makes the text, the line and the draw evolve, creating
        the frames of the animation.
        As the iteration is done for regular time, the simulation time is given by n*dt.

        Parameters
        ----------
        n : integer index for each number of iteration

        '''
        
        # update the text window 
        self.time_text.set_text('Elapsed time: {:6.2f} s'.format(n * self.dt))

        # update of each position for both the pendulums
        self.line.set_xdata([0.0, self.x1[n], self.x2[n]])
        self.line.set_ydata([0.0, self.y1[n], self.y2[n]])

        # update of each trace for both the pendulums
        if self.draw_trace:
            self.trace1.set_xdata(self.x1[0: n+1])
            self.trace1.set_ydata(self.y1[0: n+1])
            self.trace2.set_xdata(self.x2[0: n+1])
            self.trace2.set_ydata(self.y2[0: n+1])
        
        return self.line, self.trace1, self.trace2
       
    def animate(self):
        
        '''
        Function to call to generate the animation using FuncAnimation, from matplotlib.
        
        '''
        
        # cycle that iterates the sequence of frames
        self.animation = animation.FuncAnimation(fig = self.fig,
                                                 func = self.update,
                                                 frames = range(self.n_frames),
                                                 interval = 5, blit = False)