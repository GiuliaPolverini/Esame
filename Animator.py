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

    def __init__(self, pend_simulation, n_frames=-1, draw_trace = False):

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

            # tracing the second mass' phase space for θ and p
            self.theta1 = pend_simulation.integrate[:,1]
            self.p1 = pend_simulation.integrate[:,3]
            self.theta2 = pend_simulation.integrate[:,1]
            self.p2 = pend_simulation.integrate[:,3]

        if (n_frames == -1) or (n_frames > self.x1.shape[0]):
            self.n_frames = self.x1.shape[0]
        else: self.n_frames = n_frames

        l = pend_simulation.l

        self.fig, (self.ax1, self.ax2) = plt.subplots(1,2, figsize = (8,16))

        #PLOT for cartesin trajectory
        self.ax1.set_title('Cartesian Space', fontsize=10)
        #self.ax1.set_xlabel('x', fontsize=5)
        #self.ax1.set_ylabel('y', fontsize=5)
        self.ax1.set_xlim(left=-2.1*l, right=2.1*l)
        self.ax1.set_ylim(bottom=-2.1*l, top=2.1*l)

        #PLOT for phase space trajectory
        self.ax2.set_title('Phase Space', fontsize=10)

        theta1_min = np.min(self.theta1)
        if theta1_min > -np.pi: theta1_min = -np.pi
        else: theta1_min = np.floor(theta1_min/np.pi)*np.pi

        theta1_max = np.max(self.theta1)
        if theta1_max < np.pi: theta1_max = np.pi
        else: theta1_max = np.ceil(theta1_max/np.pi)*np.pi

        p1_min = np.min(self.p1)
        p1_max = np.max(self.p1)
        
        if np.isclose(p1_min, p1_max, rtol=1e-05, atol=1e-08):
            p1_min = p1_min - np.abs(p1_min)/2
            p1_max = p1_max + np.abs(p1_max)/2
        self.ax2.set_ylim(bottom=p1_min, top=p1_max)
        
        self.ax2.set_xlim(left=theta1_min,  right=theta1_max) # theta2
        self.ax2.set_ylim(bottom=p1_min*1.1, top=p1_max*1.1)
        
        theta2_min = np.min(self.theta2)
        if theta2_min > -np.pi: theta2_min = -np.pi
        else: theta2_min = np.floor(theta2_min/np.pi)*np.pi

        theta2_max = np.max(self.theta2)
        if theta2_max < np.pi: theta2_max = np.pi
        else: theta2_max = np.ceil(theta2_max/np.pi)*np.pi

        p2_min = np.min(self.p2)
        p2_max = np.max(self.p2)

        if np.isclose(p2_min, p2_max, rtol=1e-05, atol=1e-08):
            p2_min = p2_min - np.abs(p2_min)/2
            p2_max = p2_max + np.abs(p2_max)/2
        self.ax2.set_ylim(bottom=p2_min, top=p2_max)

        self.ax2.set_xlim(left=theta2_min,  right=theta2_max) # theta2
        self.ax2.set_ylim(bottom=p2_min*1.1, top=p2_max*1.1)
        #self.ax2.set_ylim(bottom=-20, top=20) # p2

        # prepare a text window for the timer
        self.time_text = self.ax1.text(0.05, 0.95,'Elapsed time: {:6.2f} s'.format(0.0),
            		                   horizontalalignment = 'left',
            		                   verticalalignment = 'top',
            		                   transform = self.ax1.transAxes)

        # the body of the pendulum is initialized by plotting the first
        # positions of the masses
        self.line, = self.ax1.plot([0.0, self.x1[0], self.x2[0]],
                                   [0.0, self.y1[0], self.y2[0]],
                                   marker = 'o')

        self.dot, = self.ax2.plot([self.theta1[0], self.p1[0]],
                                  [self.theta2[0], self.p2[0]],
                                   marker = 'o')

        # trace the whole trajectory of both the masses
        if self.draw_trace:
            self.trace1, = self.ax1.plot(self.x1[0], self.y1[0])
            self.trace2, = self.ax1.plot(self.x2[0], self.y2[0])
            self.phase_trace1, = self.ax2.plot(self.theta1[0], self.p1[0])
            self.phase_trace2, = self.ax2.plot(self.theta2[0], self.p2[0])

        else:
            self.trace1 = None
            self.trace2 = None
            self.phase_trace1 = None
            self.phase_trace2 = None


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

        self.dot.set_xdata([self.theta1[n], self.theta2[n]])
        self.dot.set_ydata([self.p1[n], self.p2[n]])

        # update of each trace for both the pendulums
        if self.draw_trace:
            self.trace1.set_xdata(self.x1[0: n+1])
            self.trace1.set_ydata(self.y1[0: n+1])

            self.trace2.set_xdata(self.x2[0: n+1])
            self.trace2.set_ydata(self.y2[0: n+1])

            self.phase_trace1.set_xdata(self.theta1[0: n+1])
            self.phase_trace1.set_ydata(self.p1[0: n+1])
            
            self.phase_trace2.set_xdata(self.theta2[0: n+1])
            self.phase_trace2.set_ydata(self.p2[0: n+1])

        return self.line, self.dot, self.trace1, self.trace2,
        self.phase_trace1, self.phase_trace2

    def animate(self):

        '''
        Function to call to generate the animation using FuncAnimation, from matplotlib.

        '''

        # cycle that iterates the sequence of frames
        self.animation = animation.FuncAnimation(fig = self.fig,
                                                 func = self.update,
                                                 frames = range(self.n_frames),
                                                 interval = 5, blit = False)
