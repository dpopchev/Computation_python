#!/usr/bin/env python3

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

# Parameters
gamma = 1
beta = 5
n = 2
args = (beta, gamma, n)

# Initial condition
ab0 = np.array([1, 1.1])

def toggle(ab, t, beta, gamma, n):
    """Right hand side for toggle ODEs."""
    a, b = ab
    return np.array([beta / (1 + b**n) - a,
                     gamma * (beta / (1 + a**n) - b)])

def plot_flow_field(ax, f, u_range, v_range, args=(), n_grid=40):
    """
    Plots the flow field with line thickness proportional to speed.

    Parameters
    ----------
    ax : Matplotlib Axis instance
        Axis on which to make the plot
    f : function for form f(y, t, *args)
        The right-hand-side of the dynamical system.
        Must return a 2-array.
    u_range : array_like, shape (2,)
        Range of values for u-axis.
    v_range : array_like, shape (2,)
        Range of values for v-axis.
    args : tuple, default ()
        Additional arguments to be passed to f
    n_grid : int, default 100
        Number of grid points to use in computing
        derivatives on phase portrait.

    Returns
    -------
    output : Matplotlib Axis instance
        Axis with streamplot included.
    """

    # Set up u,v space
    u = np.linspace(u_range[0], u_range[1], n_grid)
    v = np.linspace(v_range[0], v_range[1], n_grid)
    uu, vv = np.meshgrid(u, v)

    # Compute derivatives
    #u_vel = np.empty_like(uu)
    #v_vel = np.empty_like(vv)
    #for i in range(uu.shape[0]):
    #    for j in range(uu.shape[1]):
    #        u_vel[i,j], v_vel[i,j] = f(np.array([uu[i,j], vv[i,j]]), None, *args)

    vf = np.vectorize(f)
    u_vel, v_vel = f((uu, vv), None, *args)

    # Compute speed
    speed = np.sqrt(u_vel**2 + v_vel**2)

    # Make linewidths proportional to speed,
    # with minimal line width of 0.5 and max of 3
    lw = 0.5 + 2.5 * speed / speed.max()

    # Make stream plot
    #ax.streamplot(uu, vv, u_vel, v_vel, linewidth=lw, arrowsize=1.2,
    #              density=1, color='thistle')

    ax.quiver(u_vel, v_vel)

    return ax

if __name__ == "__main__":
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    ax.set_aspect('equal')

    ax = plot_flow_field(ax, toggle, (0, 6), (0, 6), args=args)

    plt.show()
