import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt
import matplotlib.animation as animation

#%%
def Pendulum_model(state, time, l):
    #TO DO: add l1, l2, m1, m2 as parametrs
    a1, a2, p1, p2 = state #a -> angle, p -> momentum
    g = 9.81

    expr1 = np.cos(a1 - a2)
    expr2 = np.sin(a1 - a2)
    expr3 = (1 + expr2**2)
    expr4 = p1 * p2 * expr2 / expr3
    expr5 = (p1**2 + 2 * p2**2 - p1 * p2 * expr1) * np.sin(2 * (a1 - a2)) / (2 * expr3**2)
    expr6 = expr4 - expr5

    δa1 = (p1 - p2 * expr1) / expr3
    δa2 = (2 * p2 - p1 * expr1) / expr3
    δp1 = (-2 * g * l * np.sin(a1) - expr6)
    δp2 = (-1 * g * l * np.sin(a2) + expr6)

    return δa1, δa2, δp1, δp2

#%%-----------------------------------------------------------------------------
time = np.linspace(0, 10, 2**11+1) # 2**13 = 8193 ~ 10k insanti della simulazione precedente
#        (a1,    a2,      p1,  p2)
state0 = (np.pi, np.pi-0.001, 0.0, 0.0)
l = 1.0
res = odeint(Pendulum_model, y0=state0, t=time, args=(l,))

#%%-----------------------------------------------------------------------------

a1_hat, a2_hat, p1_hat, p2_hat = res.T

x1 =  l * np.sin(a1_hat)
y1 = -l * np.cos(a1_hat)
x2 = x1 + l * np.sin(a2_hat)
y2 = y1 - l * np.cos(a2_hat)
#%%-----------------------------------------------------------------------------
#PLOTTING
# fig = plt.figure(figsize=(15,10))
# plt.plot(x2, y2, label='Pendulum 2')
# plt.legend(fontsize=20)
#%%-----------------------------------------------------------------------------

class Animator:
    def __init__(self, x1_traj, y1_traj, x2_traj, y2_traj, draw_trace = False):
        self.draw_trace = draw_trace
        self.time = 0.0

        self.x1 = x1_traj
        self.x2 = x2_traj
        self.y1 = y1_traj
        self.y2 = y2_traj

        # set up the figure
        self.fig, self.ax = plt.subplots( figsize=(12,8))
        self.ax.set_ylim(-2.5, 2.5) # variare!
        self.ax.set_xlim(-2.5, 2.5)

        # prepare a text window for the timer
        self.time_text = self.ax.text(0.05, 0.95, '',
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

    def update(self, n):
        #self.time_text.set_text('Elapsed time: {:6.2f} s'.format(self.time))

        self.line.set_xdata([0.0, self.x1[n], self.x2[n]])
        self.line.set_ydata([0.0, self.y1[n], self.y2[n]])

        if self.draw_trace:
            self.trace.set_xdata(self.x2[0: n+1])
            self.trace.set_ydata(self.y2[0: n+1])
        return self.line, self.trace

    def animate(self):
        self.animation = animation.FuncAnimation(self.fig, self.update, frames=100, interval=25, blit=False)

#%%
animator = Animator(x1, y1, x2, y2, draw_trace = True)
animator.animate()
plt.show()
