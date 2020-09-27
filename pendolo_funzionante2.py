# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 14:01:55 2020

@author: Giulia
"""

import matplotlib
matplotlib.use('Qt4Agg') #'tkAgg' if Qt not present
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
class Pendulum:
    def __init__(self, a1, a2, dt, length = 1.0):
        self.a1 = a1
        self.a2 = a2
        self.l = length
        self.p1 = 0.0
        self.p2 = 0.0
        self.dt = dt
        self.g = 9.81
        first_position = self.polar_to_cartesian()
        self.trajectory = [first_position]
    def polar_to_cartesian(self):
        x1 =   self.l * np.sin(self.a1)
        y1 = - self.l * np.cos(self.a1)
        x2 = x1 + self.l * np.sin(self.a2)
        y2 = y1 - self.l * np.cos(self.a2)
        return np.array([[0.0,0.0], [x1,y1], [x2,y2]])
    def evolve(self):
        a1 = self.a1
        a2 = self.a2
        p1 = self.p1
        p2 = self.p2
        g = self.g
        l = self.l
        eq1 = np.cos(a1 - a2)
        eq2 = np.sin(a1 - a2)
        eq3 = (1 + eq2**2)
        eq4 = p1 * p2 * eq2 / eq3
        eq5 = (p1**2 + 2 * p2**2 - p1 * p2 * eq1) \
                *np.sin(2 * (a1 - a2)) / 2 / eq3**2
        eq6 = eq4 - eq5
        a1 += self.dt * (p1 - p2 * eq1) / eq3
        a2 += self.dt * (2 * p2 - p1 * eq1) / eq3
        p1 += self.dt * (-2 * g * l * np.sin(a1) - eq6)
        p2 += self.dt * (    -g * l * np.sin(a2) + eq6)
        self.a1 = a1
        self.a2 = a2
        self.p1 = p1
        self.p2 = p2
        new_position = self.polar_to_cartesian()
        self.trajectory.append(new_position)
        return new_position
class Animator:
    def __init__(self, pendulum, draw_trace = False):
        self.pendulum = pendulum
        self.d_t = draw_trace
        self.time = 0.0
        #set up the figure
        self.fig, self.ax = plt.subplots()
        self.ax.set_ylim(-2.5, 2.5)
        self.ax.set_xlim(-2.5, 2.5)
        #prepare a text window for the timer
        self.time_text = self.ax.text(0.05, 0.95, '', horizontalalignment = 'left', verticalalignment = 'top', transform = self.ax.transAxes)
        #initialize by plotting the last position of the trajectory
        self.line, = self.ax.plot(self.pendulum.trajectory[-1][:, 0], self.pendulum.trajectory[-1][:, 1], marker = 'o')
        #trace the whole trajectory of the second mass
        if self.d_t:
            self.trace, = self.ax.plot([a[2, 0] for a in self.pendulum.trajectory], [a[2, 1] for a in self.pendulum.trajectory])
    def advance_time_step(self):
        while True:
            self.time += self.pendulum.dt
            yield self.pendulum.evolve()
    def update(self, data):
        self.time_text.set_text('Time:{:6.2f}s'.format(self.time))
        self.line.set_ydata(data[:, 1])
        self.line.set_xdata(data[:, 0])
        if self.d_t:
            self.trace.set_xdata([a[2,0] for a in self.pendulum.trajectory])
            self.trace.set_ydata([a[2,1] for a in self.pendulum.trajectory])
        return self.line,
    def animate(self):
        self.animation = animation.FuncAnimation(self.fig, self.update, self.advance_time_step, interval=25, blit=False)
pendulum = Pendulum(a1 = np.pi, a2 = np.pi/2 - 0.01, dt = 0.01)







# animator = Animator(pendulum = pendulum, draw_trace = True)
# animator.animate()
# plt.show()
