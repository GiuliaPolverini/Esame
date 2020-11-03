# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:56:19 2020

@author: Giulia
"""
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation

class Animator: # tutto ciò che concerne l'animazione
                # per funzionare ha bisogno della simulazione
    
    def __init__(self, pend_simulation, draw_trace = False):
        self.draw_trace = draw_trace
        self.dt = pend_simulation.dt # accedo all'attributo dt dell'oggetto 
                                     # pend_simulation 
        # dt è la durata in sec che c'è tra un tempo e l'altro della simulazione
        # (lo calcolo nella classe PendulumSimulation facendo la differenza tra
        # 2 istanti successivi)

        if pend_simulation.motion is None:
            print('Moto del pendolo non simulato.')
        else:
            # Saving the trajectories (queste sono le serie di posiz delle 2m)
            self.x1 = pend_simulation.x1
            self.y1 = pend_simulation.y1
            self.x2 = pend_simulation.x2
            self.y2 = pend_simulation.y2 # occupo questi valori creando un nuovo 
                                         # attributo stavolta della class Animator

            self.bound = 0.2 #boundary percentage: serve per il plot, insieme 
            # al pezzo sotto adatta lo spazio della rappresentazione alla 
            # lunghezza + il 20%

            x2_extr = np.asarray([np.min(self.x2), np.max(self.x2)]) * (1 + self.bound)
            y2_extr = np.asarray([np.min(self.y2), np.max(self.y2)]) * (1 + self.bound)
            # queste sono le coord degli estremi da dare alla figure; ho trasformato 
            # la lista tra le [] dove ci sono i min e i max di x e y in un numpy 
            # array, perché così posso moltiplicare tutto l'array per un singolo 
            # numero (1+0.2 qui, dunque il 120%) --> numpy capisce che deve 
            # moltiplicare per ogni singolo elemento

            # set up the figure
            self.fig, self.ax = plt.subplots(figsize = (8,8))
            self.ax.set_xlim(*x2_extr) #self.ax.set_xlim(x2_extr[0], x2_extr[1])
            self.ax.set_ylim(*y2_extr) 
            # l'asterisco passa come singoli parametri alle funzioni set_xlim
            # e set_ylim tutti i valori contenuti nella sequenza --> ciascuno
            # passa quindi 2 param in questo caso, [0] e [1]

            # prepare a text window for the timer
            self.time_text = self.ax.text(0.05, 0.95, 'Elapsed time: {:6.2f} s'.format(0.0),
                		                  horizontalalignment = 'left',
                		                  verticalalignment = 'top',
                		                  transform = self.ax.transAxes)

            # initialize by plotting the last position of the trajectory
            self.line, = self.ax.plot([0.0, self.x1[0], self.x2[0]],
                                      [0.0, self.y1[0], self.y2[0]],
                                      marker = 'o')
            # la line è il corpo del pendolo; poi serie di valori per x e y da
            # accoppiare in 3 coppie che danno i 3 punti iniziali con marker
            # tondo, che danno la linea spezzata

            # trace the whole trajectory of both the masses, inizialmente
            # col primo punto della traiettoria
            if self.draw_trace:
                self.trace1, = self.ax.plot(self.x1[0], self.y1[0])
                self.trace2, = self.ax.plot(self.x2[0], self.y2[0])
            else:
                self.trace = None

    def update(self, n):
        self.time_text.set_text('Elapsed time: {:6.2f} s'.format(n * self.dt))
        # cambia il testo: il format prende quello che sta nelle tonde e lo
        # mette nelle graffe in cui dico che voglio al max 6 cifre davanti alla
        # virgola, al max 2 dopo, ed è un numero float
        # Do n perché lui sa che sto facendo l'n-esima iterazione (lo sa da 
        #range) e visto che l'iterazione è fatta a tempi regolari allora il 
        # tempo simulato è n*dt

        # self.time_text.set_text('Frame: {} '.format(n))
        self.line.set_xdata([0.0, self.x1[n], self.x2[n]])
        self.line.set_ydata([0.0, self.y1[n], self.y2[n]])
        # modifica la posizione del pendolo, perché sovrascrive le coord x, y

        if self.draw_trace:
            self.trace1.set_xdata(self.x1[0: n+1])
            self.trace1.set_ydata(self.y1[0: n+1])
            self.trace2.set_xdata(self.x2[0: n+1])
            self.trace2.set_ydata(self.y2[0: n+1])
        return self.line, self.trace1, self.trace2
        # sovrascrive la traccia, prende i primi n valori della traiettoria
        # n+1 per lo slicing di numpy ??

    def animate(self):
        self.animation = animation.FuncAnimation(fig = self.fig,
                                                 func = self.update,
                                                 frames = range(1000),
                                                 interval = 5, blit = False)
        # metodo animate, che con FuncAnimation genera un frame per ogni
        # valore che c'è nella sequenza frames, prende ognuno di questi valori
        # e li passa a func, che è la responsabile dell'aggiornamento del nuovo
        # frame nella fig scelta.
        # Dunque è un ciclo che itera la sequenza di frames, perciò undate
        # viene eseguita 1000 volte qui e genera una sequenza di 1000 oggetti