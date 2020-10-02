import matplotlib.pyplot as plt
import numpy as np

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

#%%

pendulum = Pendulum(theta1=np.pi, theta2=np.pi/2, dt=0.01)

n_evolution = 999

#sviluppo moto singolo pendolo
for i in range(0, n_evolution):
    pendulum.evolve()

#conversione di una lista di np.array in 2d in un unico np.array a 3d
traj_tensor = np.asarray(pendulum.trajectory)

#1 dim: istante della traiettoria
#2 dim: 0: perno, 1:massa di mezzo, 2:estremità
#3 dim: 0: x, 1:y

#%%
#sottocampionamento della traiettoria prendendo i dispari
odd_index = range(1, 999, 2) # [1, 3, 5,..., 999]
sub_sampled_traj = traj_tensor[odd_index,:,:]

#%%
#Codice per plottare gli andamenti delle 2 masse (insieme)
fig = plt.figure()

for mass in [1, 2]:
    plt.plot(trajectory[:,mass,0], 
             trajectory[:,mass,1])
    
#plt.savefig('Esame/2masses.png')
#%%
pendulum2 = Pendulum(theta1=np.pi, theta2=np.pi/2 - 0.01, dt=0.01)
pendulum3 = Pendulum(theta1=np.pi, theta2=np.pi/2 - 0.02, dt=0.01)
for i in range(0, n_evolution):
    pendulum2.evolve()
    pendulum3.evolve()
    
trajectory2 = np.asarray(pendulum2.trajectory)
trajectory3 = np.asarray(pendulum3.trajectory)
mass = 2


plt.plot(trajectory[:, mass, 0],
         trajectory[:, mass, 1])
plt.plot(trajectory2[:, mass, 0],
         trajectory2[:, mass, 1])
plt.plot(trajectory3[:, mass, 0],
         trajectory3[:, mass, 1])
#%%

delta_t2 = np.linspace(0, 0.99, 100)
traj_list = []
mass = 2

for d in delta_t2:
    pend = Pendulum(theta1 = np.pi, theta2 = np.pi/2 - d, dt = 0.01)
    
    #per ognuno dei 100 pendoli "creati sopra" si fanno 1000 evoluzioni
    for i in range(0, n_evolution):
        pend.evolve()
    
    traj = np.asarray(pend.trajectory) #conversione
    traj_mass2 = traj[:, mass, :] #slicing 2°massa
    traj_list.append(traj_mass2) #salvataggio

for traj in traj_list:
    plt.plot(traj[:,0],traj[:,1])
    
#%%









