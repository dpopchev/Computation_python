#!/usr/bin/env python3
# do not hesitate to debug
import pdb

# python computation modules and visualization
import numpy as np
import sympy as sy
import scipy as sp
import matplotlib.pyplot as plt

from sympy import Q as syQ

sy.init_printing(use_latex=True,forecolor="White")

def Lyapunov_stability_test_linear(ev):
    ''' test if a linear homogeneous system with constant coefficients is stable
    in the sense of Lyapunov by checking the theorem conditions against the
    provided eigenvalues
    source https://www.math24.net/stability-theory-basic-concepts/
    TODO taking into account eigenvalue multiplicity '''

    # the criteria result will be saved here
    r = None

    # system is asymptotically stable if only if
    # all eigenvalues have negative real parts
    r = 'asymptotically stable' if ( not r
            and  all(sy.ask(syQ.negative(sy.re(_))) for _ in ev) ) else None

    # system is stable if and only if
    # all eigenvalues have nonpositive real parts
    # TODO incorporate algebraic and geometric multiplicity criteria
    r = 'stable' if ( not r
            and all(sy.ask(syQ.nonpositive(sy.re(_))) for _ in ev) ) else None

    # system is unstable if
    # at least one eigenvalue has positive real part
    # TODO incorporate algebraic and geometric multiplicity criteria
    r = 'unstable' if ( not r
            and any(sy.ask(syQ.positive(sy.re(_))) for _ in ev) ) else None

    return r

def Lyapunov_stability_test_nonlinear(ev):
    ''' test if the fixed point of a nonlinear structure stable system
    is stable, unstable, critical or impossible to determine using Lyapunov
    criteria of first order and thus other methods are needed
    TODO tests are only applicable for structurally stable systems, i.e.
    with purely imaginary eigenvalues are not taken into account
    source https://www.math24.net/stability-first-approximation/ '''

    # the criteria result will be saved here
    r = None

    # system is asymptotically stable if only if
    # all eigenvalues have negative real parts
    r = 'asymptotically stable' if ( not r
            and  all(sy.ask(syQ.negative(sy.re(_))) for _ in ev) ) else None

    # system is unstable if
    # at least one eigenvalue has positive real part
    r = 'unstable' if ( not r
            and any(sy.ask(syQ.positive(sy.re(_))) for _ in ev) ) else None

    # if all eigenvalues have non-positive real parts,
    # and there is at least one eigenvalue with zero real part
    # then fixed point can be stable or unstable and other methods should be
    # used, thus mark the point critical
    r = 'critical' if ( not r
            and all(sy.ask(Q.nonpositive(sy.re(_))) for _ in ev)
            and any(sy.re(_) == 0 for _ in ev)
            ) else None

    return r if r else 'not decided'

def RouthHurwitz_Criterion(p):
    ''' return principal minors of Hurwitz matrix as sympy polynomials, which if
    all are positive it is sufficient condition for asymptotic stability
    NOTE: if all n-1 principal minors are positive, and nth minor is zero,
    the system is at the boundary of stability, with two cases:
        a_n = 0 -- one of the root is zero and system is on the boundary of
        aperiodic stability
        n-1 minor is zero -- there are two complex conjugate imaginary roots and
        the system is at boundary of oscillatory stability
    source https://www.math24.net/routh-hurwitz-criterion/ '''

    # initial key and index pair needed to create Hurwitz matrix via sympy banded
    # each entry is of the type [ dictionary key, coefficient slice ]
    idxs = [ [ 1, 0 ] ]

    # generate next key by decrementing with 1
    genKey = lambda _: _ - 1

    # generate next index by incrementing with 1 if key was nonnegative
    # or with 2 if key is negative
    genSlice = lambda _, __: __ + 1 if _ >= 0 else __ + 2

    # fill the rest pairs w.r.t. the polynomial degree - 1, as we already have
    # one entry
    for _ in range(p.degree() - 1):
        key = genKey(idxs[-1][0])
        idxs.append( [ key, genSlice(key, idxs[-1][1] ) ] )

    # create the matrix itself
    H = sy.banded({ k: p.all_coeffs()[v:] for k, v in idxs })

    return [ H[:_, :_].det() if _ > 0 else p.LC() for _ in range(0, p.degree()+1) ]

# define independent variable
t = sy.symbols('t', real=True)

# define dependent variables individually and pact them in an variable
theta, omega = sy.symbols(r'\theta, \omega', real = True)
Y = theta, omega

# define free parameters of they system and pack them in a variable
g, L =  sy.symbols('g, L', positive = True)
parms = g, L

# create rhs as sympy expressions
theta_dt = omega
omega_dt = -(g/L)*sy.sin(theta)
rhs = {}
rhs['sympy'] = sy.Matrix([theta_dt, omega_dt])

# convert the sympy matrix function to numpy function with usual signature
rhs['numpy'] = sy.lambdify((t, Y, *parms), rhs['sympy'], 'numpy')

# create Jacobian matrix as sympy expression
J = {}
J['sympy'] = rhs['sympy'].jacobian(Y)

# convert the sympy Jacobian expression to numpy function with usual signature
J['numpy'] = sy.lambdify((t, Y, *parms), J['sympy'])

# calculate rhs fixed points
fixed_points = sy.solve(rhs['sympy'], Y)

# substitute each fixed point in the Jacobian
# and calculate the eigenvalues
J_fixed = {}
for i, fp in enumerate(fixed_points):

    J_subs = J['sympy'].subs( [(y, v) for y, v in zip(Y, fp)])

    #J_eigenvals = J_subs.eigenvals(multiple=True)
    J_eigenvals = J_subs.eigenvals()

    # save the fixed point results in more details
    # most importantly the eigenvalues and their corresponding multiplicity
    J_fixed[i] = {
            'fixed point': fp,
            'subs': J_subs,
            'eigenvalues': list(J_eigenvals.keys()),
            'multiplicity': list(J_eigenvals.values())
            }

def plot_phase_portrait(ax, rhs, section, args=(), n_points=25):
    ''' plot section of phase space of a field defined via its rhs '''

    # create section grid
    x_grid, y_grid = np.meshgrid(
            np.linspace( section[0][0], section[0][1], n_points ),
            np.linspace( section[1][0], section[1][1], n_points )
            )

    # calculate rhs on the grid
    xx, yy = rhs(None, ( x_grid, y_grid ), *args)

    # compute vector norms and make line width proportional to them
    # i.e. greater the vector length, the thicker the line
    # TODO not sure why rhs returns different shape
    vector_norms = np.sqrt(xx[0]**2 + yy[0]**2)
    lw = 0.25 + 3*vector_norms/vector_norms.max()

    # plot the phase portrait
    ax.streamplot(
            x_grid, y_grid,
            xx[0], yy[0],
            linewidth = lw,
            arrowsize = 1.2,
            density = 1
            )

    return ax

def plot_main():

    fig, ax = plt.subplots()

    ax = plot_phase_portrait(
            ax,
            rhs['numpy'],
            (
                ( -np.pi, np.pi ),
                ( -2*np.pi, 2*np.pi)
                ),
            args = ( 5, 1 ),
            )

if __name__ == '__main__':

    plot_main()
