# Esame
### The Double Pendulum
![nome_immagine](https://physicspython.files.wordpress.com/2019/02/double_pendulum-e1549896953882.png?w=600)

The double pendulum is a classical exercise of mechanics. Consider the approximation that each pendulum consists of a point mass ![equation1](https://latex.codecogs.com/gif.latex?m) hanging on an ideal (non-elastic, mass-less) string of length ![equation2](https://latex.codecogs.com/gif.latex?l) in a constant, homogeneous  gravitational field of acceleration ![equation3](https://latex.codecogs.com/gif.latex?g). While the first pendulum is attached to a rigid, motionless point, the second is attached to the point mass of the first one,  drawing a chaotic trajectory. Moreover, the motion here only happens in a single plane, so that each pendulum is attached by a rigid string to a specific point, and to greatly semplify the equation I assumed that ![equation4](https://latex.codecogs.com/gif.latex?m_{1}&space;=&space;m_{2}&space;=&space;1) and ![equation5](https://latex.codecogs.com/gif.latex?l_{1}&space;=&space;l_{2}&space;=&space;1.)

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

1) First the user has to choose the configurations he/she prefers, using the syntax of [config](https://github.com/GiuliaPolverini/Esame/blob/master/config.json);
2) Then to make the simulation run the user has to FILE DOUBLEPEND
3) At the end ...

This is how I divided my project into blocks:

- In the file [PendulumSimulation](https://github.com/GiuliaPolverini/Esame/blob/master/PendulumSimulation.py) I have initialized some parameters as the initial state of the system and the time, solved the differential equations for the system thanks to the odeint package and got the cartesian coordinates for each position; in the last part I tested ....

- In the file [Animator](https://github.com/GiuliaPolverini/Esame/blob/master/Animator.py) I have set the length of time for the simulation, the size of the animation window and of the plot and how the text window has to evolve; then I have built the line of the pendulum with its three points (the central pin and the two masses), inizializing their positions and the trace to draw the trajectory; at the end I created a cycle to iterate the frames and update the positions and the traces.

- In the file [config](https://github.com/GiuliaPolverini/Esame/blob/master/config.json) there are the parameters the user can vary according to how he/she wants them.

- In the file [DoublePendulum](https://github.com/GiuliaPolverini/Esame/blob/master/DoublePendulum.py) first I have brought out the values I needed from the [config](https://github.com/GiuliaPolverini/Esame/blob/master/config.json) file to create the system; then I have given the comand to create the simulation and the animation.

This is how the animation looks like after 9,72 seconds:
![ScreenShot](Screenshot.png)
