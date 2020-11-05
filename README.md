# Esame
### The Double Pendulum
![nome_immagine](https://physicspython.files.wordpress.com/2019/02/double_pendulum-e1549896953882.png?w=600)

The double pendulum is a classical exercise of classical mechanics. I made the approximation that each pendulum consists of a point mass ![equation1](https://latex.codecogs.com/gif.latex?m) hanging on an ideal (non-elastic, mass-less) string of length ![equation2](https://latex.codecogs.com/gif.latex?l) in a constant, homogeneous  gravitational field of acceleration ![equation3](https://latex.codecogs.com/gif.latex?g). While the first pendulum is attached to a rigid, motionless point, the second is attached to the point mass of the first one,  drawing a chaotic trajectory. Moreover, I assumed that the motion only happens in a single plane, so that each pendulum is attached by a rigid string to a specific point, and to greatly semplify the equation I set ![equation4](https://latex.codecogs.com/gif.latex?m_{1}&space;=&space;m_{2}&space;=&space;1) and ![equation5](https://latex.codecogs.com/gif.latex?l_{1}&space;=&space;l_{2}&space;=&space;1.)
This means that the system is described by a set of four first order differential equations, two for each degree of freedom:

![equation6](https://latex.codecogs.com/gif.latex?\dot{\theta_1}&space;=&space;\frac{p_1&space;-&space;p_2cos(\theta_1&space;-&space;\theta_2)}{1&space;+&space;sin^2(\theta_1&space;-&space;\theta_2)})

![equation7](https://latex.codecogs.com/gif.latex?\dot{\theta_2}&space;=&space;\frac{2p_2&space;-&space;p_1cos(\theta_1&space;-&space;\theta_2)}{1&space;+&space;sin^2(\theta_1&space;-&space;\theta_2)})

![equation8](https://latex.codecogs.com/gif.latex?\dot{p_1}&space;=&space;-2gsin(\theta_1)&space;-&space;A&space;+&space;B)

![equation9](https://latex.codecogs.com/gif.latex?\dot{p_2}&space;=&space;-gsin(\theta_2)&space;+&space;A&space;-&space;B)

where

![equation10](https://latex.codecogs.com/gif.latex?A&space;=&space;\frac{p_1p_2sin(\theta_1&space;-&space;\theta_2)}{1&space;+&space;sin^2(\theta_1&space;-&space;\theta_2)})

![equation11](https://latex.codecogs.com/gif.latex?B&space;=&space;\frac{p_1^2&space;+&space;2p_2^2&space;-&space;p_1p_2cos(\theta_1&space;-&space;\theta_2)}{2[1&space;+&space;sin^2(\theta_1&space;-&space;\theta_2)]^2}sin[2(\theta_1&space;-&space;\theta_2)])


### Structure of the project
These are the steps in order to start the program and get the animation:

1) The user has to choose the configurations he/she prefers, using the syntax of [config](https://github.com/GiuliaPolverini/Esame/blob/master/config.json); in particular, the user has to specify the two angles and the two angular momenta for both the systems ...
2) ...
3) ...

This is how I divided my project into blocks:

- In the file [PendulumSimulation](https://github.com/GiuliaPolverini/Esame/blob/master/PendulumSimulation.py) I have built the ...

- In the file [Animator](https://github.com/GiuliaPolverini/Esame/blob/master/Animator.py) ...

- In the fine [DoublePendulum](https://github.com/GiuliaPolverini/Esame/blob/master/DoublePendulum.py) ...

- In the file [config](https://github.com/GiuliaPolverini/Esame/blob/master/config.json) there are the parameters the user can vary according to COME LA VUOLE VEDERE!


----------------------------------------------------------------------------------------------------------------------------------------------------
{ising: I have built the Ising functions that calculate energy and magnetization of the system, and save them in arrays in order to use them for further data analysis. In addition, for a given temperature there is a function that calculates the different states of a lattice during time, from disordered state to ordered state.}

{testing: I have tested all the Ising functions to ensure that all of them work properly, using hypothesis testing.}

{configuration: there are all the definitions of the parameters used in the simulation file, as number of spins per lattice (N*M), temperature intervals and so on. Furthermore, there are the local paths in order to load the array data and to save them as images and graphs. It's a .txt file that is imported in simulation file.}

{simulation: there is the main part of the code, where I have used the functions of ising file in order to calculate the energy and the magnetization of a configuration of spins for a range of temperatures across the critical one ***T<sub>c***, showing a steeply decrease in energy from high temperatures to low ones and a rapidly increase in magnetization, a clear sign of a phase transition. In addition there is the calculation of the different states of the configuration of spins for a given temperatrue, lower than ***T<sub>c***, respect to time, which shows that the system coarsens toward the configuration of all spins aligned; then I saved these states in an array to process them in further data analysis.}
