# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 18:30:29 2020

@author: Giulia
"""

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
        p2 += self.dt * (-g * l * np.sin(theta2) + expr6)

        self.theta1 = theta1
        self.theta2 = theta2
        self.p1 = p1
        self.p2 = p2

        new_position = self.polar_to_cartesian()
        self.trajectory.append(new_position)
        return new_position