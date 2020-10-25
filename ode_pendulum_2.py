import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt
import matplotlib.animation as animation

class Animator:
    # TO ADD: 'dt' as parameter
    def __init__(self, pend_simulation, draw_trace = False):
        self.draw_trace = draw_trace
        self.dt = pend_simulation.dt

        # Saving the trajectories
        self.x1 = pend_simulation.x1
        self.y1 = pend_simulation.y1
        self.x2 = pend_simulation.x2
        self.y2 = pend_simulation.y2

        self.bound = 0.2 #boundary percentage

        x2_extr = np.asarray([np.min(self.x2), np.max(self.x2)]) * (1 + self.bound)
        y2_extr = np.asarray([np.min(self.y2), np.max(self.y2)]) * (1 + self.bound)

        # set up the figure
        self.fig, self.ax = plt.subplots( figsize=(12,8))
        self.ax.set_xlim(*x2_extr)
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

class PendulumSimulation:
    def __init__(self, state0, time_end, exp_time_sampling = 11):
        self.time = np.linspace(0.0, time_end, 2**exp_time_sampling+1)
        self.dt = self.time[1] - self.time[0]
        self.l = 1
        self.state0 = state0
        self.motion = None

    def Pendulum_model(self, state, time, l):
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

    def simulate(self):
        self.motion = odeint(self.Pendulum_model, y0=self.state0, t=self.time, args=(self.l,))
        self.cartesian_traj()

    def cartesian_traj(self):
        a1_hat, a2_hat, p1_hat, p2_hat = self.motion.T
        self.x1 =  self.l * np.sin(a1_hat)
        self.y1 = -self.l * np.cos(a1_hat)
        self.x2 = self.x1 + self.l * np.sin(a2_hat)
        self.y2 = self.y1 - self.l * np.cos(a2_hat)

pend_sim = PendulumSimulation(time_end = 40, state0 = (np.pi, np.pi/2, 0.0, 5.0))
pend_sim.simulate()
anim = Animator(pend_sim, draw_trace = True)
anim.animate()
plt.show()

#%%-----------------------------------------------------------------------------
# DERIVATIVE FUNCTION FOR THE DOUBLE PENDULUM SYSTEM -> attribute for OdeInt
# def Pendulum_model(state, time, l):
#     #TO DO: add l1, l2, m1, m2 as parametrs
#     a1, a2, p1, p2 = state #a -> angle, p -> momentum
#     g = 9.81
#
#     expr1 = np.cos(a1 - a2)
#     expr2 = np.sin(a1 - a2)
#     expr3 = (1 + expr2**2)
#     expr4 = p1 * p2 * expr2 / expr3
#     expr5 = (p1**2 + 2 * p2**2 - p1 * p2 * expr1) * np.sin(2 * (a1 - a2)) / (2 * expr3**2)
#     expr6 = expr4 - expr5
#
#     δa1 = (p1 - p2 * expr1) / expr3
#     δa2 = (2 * p2 - p1 * expr1) / expr3
#     δp1 = (-2 * g * l * np.sin(a1) - expr6)
#     δp2 = (-1 * g * l * np.sin(a2) + expr6)
#
#     return δa1, δa2, δp1, δp2
#%%-----------------------------------------------------------------------------
#SETTING FOR ODEINT
# exp_time_sampling = 11
# time_end = 40
# time = np.linspace(0.0, t_end, 2**exp_time_sampling+1) # 2**13 = 8193 ~ 10k insanti della simulazione precedente
#        (a1,    a2,      p1,  p2)
# state0 = (np.pi, np.pi/2, 0.0, 5.0)
# l = 1.0
#%%-----------------------------------------------------------------------------
# SIMULATION: Call to OdeInt
# res = odeint(Pendulum_model, y0=state0, t=time, args=(l,))
#%%-----------------------------------------------------------------------------
#Conversion to cartesian coordinates
# a1_hat, a2_hat, p1_hat, p2_hat = res.T
# x1 =  l * np.sin(a1_hat)
# y1 = -l * np.cos(a1_hat)
# x2 = x1 + l * np.sin(a2_hat)
# y2 = y1 - l * np.cos(a2_hat)
# dt_pendulum = time[1] - time[0]

#%%-----------------------------------------------------------------------------
#PLOTTING
# fig = plt.figure(figsize=(15,10))
# plt.plot(x2, y2, label='Pendulum 2')
# plt.scatter(x2[0], y2[0], marker = 'o', c='r')
# plt.legend(fontsize=20)
