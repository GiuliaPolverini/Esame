import numpy as np
from scipy.integrate import odeint # package for ODE integration
import matplotlib.pylab as plt
import matplotlib.animation as animation
import json

#%%-----------------------------------------------------------------------------
class Animator:
    # TO ADD: 'dt' as parameter
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
            self.fig, self.ax = plt.subplots( figsize=(12,8))
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

class PendulumSimulation:
    '''
    Classe per la simulazione del sistema fisico del doppio pendolo...
    '''
    def __init__(self, state0, time_end, exp_time_sampling = 11, l = 1):
        '''
        Classe per la simulazione del sistema fisico del doppio pendolo...
        '''
        self.time = np.linspace(0.0, time_end, 2**exp_time_sampling+1)
        self.dt = self.time[1] - self.time[0]
        self.l = l
        self.state0 = state0
        self.motion = None

    def Pendulum_model(self, state, time, l):
        ''' (ESEMPIO DI DOCUMENTAZIONE)
        (descrizione) Funzione per la derivata del sistema....

        parametrs:
            state: stato iniziale del sistema
            time : campionamento temporale su cui fare la simulazione
                   ~ np.linspace(0,10, 2**9+1)
            l    : lunghezza delle sbarre del pendolo. Nell'approssimazione in cui l1 = l2.

        return:
            δa1, δa2, δp1, δp2: variazioni infinitesime delle 4 variabili del sistema.
        '''
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
        '''
        Function for the conversion to cartesian coordinate of the series of values in 'self.motion'.
        '''
        a1_hat, a2_hat, p1_hat, p2_hat = self.motion.T
        self.x1 =  self.l * np.sin(a1_hat)
        self.y1 = -self.l * np.cos(a1_hat)
        self.x2 = self.x1 + self.l * np.sin(a2_hat)
        self.y2 = self.y1 - self.l * np.cos(a2_hat)

#%%-----------------------------------------------------------------------------
#loading configuration
with open('./config.json') as f:
  config = json.load(f)
config['pend']['state0'] = [float(st_vl) for st_vl in config['pend']['state0']]

#pendulum simulation
pend_sim = PendulumSimulation(**config['pend']) #<- "key = value" per tutti gli elementi del dizionario
pend_sim.simulate()

#pendulum animation
anim = Animator(pend_sim, **config['anim'])
anim.animate()
plt.show()
