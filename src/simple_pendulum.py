#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def pendulum(t, y, g, L):
    """ simple pendulum ODE in polar coordinates

    parameters

        t: float
            independent parameter

        y: vector
            y[0]    angular position
            y[1]    angular velocity

        g: float
            gravitational acceleration

        L: float
            pendulum length

    returns
        dydt: vector
            r.h.s. of the ODE
    """

    thetha, omega = y

    dthethadt = omega
    domegadt = (-g/L)*np.sin(thetha)

    return np.array([ dthethadt, domegadt ])

def vector_field():
    """ quiver wrapper to create vector field of a ODE system """

    # set phase portrait boundaries
    thetha_min = -2*np.pi
    thetha_max = 2*np.pi
    thetha_step = 0.5

    omega_min = -4*np.pi
    omega_max = 4*np.pi
    omega_step = 1

    # set the free parameters of the ODE
    g, L = 9.81, 1

    # create grid where to place the vectors
    thetha, omega = np.meshgrid(
            np.arange(thetha_min, thetha_max, thetha_step),
            np.arange(omega_min, omega_max, omega_step)
            )

    # vectorize the function ODE so to take advantage of numpy WoW
    vpendulum = np.vectorize(pendulum)

    # compute the vector components and their lengths
    u, v = vpendulum(None, (thetha, omega), g, L)

    fig, ax = plt.subplots(1,1)
    ax.quiver(u, v)

    plt.show()

    return 0

if __name__ == '__main__':

    vector_field()
