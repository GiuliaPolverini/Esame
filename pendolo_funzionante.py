import matplotlib
matplotlib.use('Qt4Agg') # 'tkAgg' if Qt not present
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

class Pendulum:
    def __init__(self, theta1, theta2, dt, length = 1.0):
        self.theta1 = theta1
        self.theta2 = theta2
        self.length = length

        self.p1 = 0.0
        self.p2 = 0.0

        self.dt = dt

        self.g = 9.81

        first_position = self.polar_to_cartesian()
        self.trajectory = [first_position]

    def polar_to_cartesian(self):
        x1 =  self.length * np.sin(self.theta1)
        y1 = -self.length * np.cos(self.theta1)
        x2 = x1 + self.length * np.sin(self.theta2)
        y2 = y1 - self.length * np.cos(self.theta2)

        return np.array([[0.0, 0.0], [x1, y1], [x2, y2]])

    def evolve(self):
        theta1 = self.theta1
        theta2 = self.theta2
        p1 = self.p1
        p2 = self.p2
        g = self.g
        l = self.length

        expr1 = np.cos(theta1 - theta2)
        expr2 = np.sin(theta1 - theta2)
        expr3 = (1 + expr2**2)
        expr4 = p1 * p2 * expr2 / expr3
        expr5 = (p1**2 + 2 * p2**2 - p1 * p2 * expr1) \
        * np.sin(2 * (theta1 - theta2)) / 2 / expr3**2
        expr6 = expr4 - expr5

        theta1 += self.dt * (p1 - p2 * expr1) / expr3
        theta2 += self.dt * (2 * p2 - p1 * expr1) / expr3
        p1 += self.dt * (-2 * g * l * np.sin(theta1) - expr6)
        p2 += self.dt * (    -g * l * np.sin(theta2) + expr6)

        self.theta1 = theta1
        self.theta2 = theta2
        self.p1 = p1
        self.p2 = p2

        new_position = self.polar_to_cartesian()
        self.trajectory.append(new_position)
        return new_position

class Animator:
    def __init__(self, pendulum, draw_trace=False):
        self.pendulum = pendulum
        self.draw_trace = draw_trace
        self.time = 0.0

        # set up the figure
        self.fig, self.ax = plt.subplots()
        self.ax.set_ylim(-2.5, 2.5)
        self.ax.set_xlim(-2.5, 2.5)

        # prepare a text window for the timer
        self.time_text = self.ax.text(0.05, 0.95, '',
                                      horizontalalignment='left',
                                      verticalalignment='top',
                                      transform=self.ax.transAxes)

        # initialize by plotting the last position of the trajectory
        self.line, = self.ax.plot(self.pendulum.trajectory[-1][:, 0],
                                  self.pendulum.trajectory[-1][:, 1],
                                  marker='o')

        # trace the whole trajectory of the second pendulum mass
        if self.draw_trace:
            self.trace, = self.ax.plot([a[2, 0] for a in self.pendulum.trajectory],
                                       [a[2, 1] for a in self.pendulum.trajectory])

    def advance_time_step(self):
        while True:
            self.time += self.pendulum.dt
            yield self.pendulum.evolve()

    def update(self, data):
        self.time_text.set_text('Elapsed time: {:6.2f} s'.format(self.time))

        self.line.set_ydata(data[:, 1])
        self.line.set_xdata(data[:, 0])

        if self.draw_trace:
            self.trace.set_xdata([a[2, 0] for a in self.pendulum.trajectory])
            self.trace.set_ydata([a[2, 1] for a in self.pendulum.trajectory])
        return self.line,

    def animate(self):
        self.animation = animation.FuncAnimation(self.fig, self.update,
                         self.advance_time_step, interval=25, blit=False)

pendulum = Pendulum(theta1=np.pi, theta2=np.pi/2-0.01, dt=0.01)
len(pendulum.trajectory)

for i in range(0,999):
    pendulum.evolve()

len(pendulum.trajectory)

trajectory = np.asarray(pendulum.trajectory)
trajectory.shape

n_massa = 2
fig = plt.figure()
plt.plot(trajectory[:,n_massa,0],trajectory[:,n_massa,1])
plt.savefig('traj_2.png')
plt.show()






#1 dim: istante della traiettoria
#2 dim: 0: perno, 1:massa di mezzo, 2:estremità
#3 dim: 0: x, 1:y

#%%-----------------------------------
# animator = Animator(pendulum=pendulum, draw_trace=True)
# animator.animate()
# plt.show()
