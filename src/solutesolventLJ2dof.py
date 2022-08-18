#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 7 19:28:35 2020

@author: Shibabrat Naik, shiba@vt.edu
"""

import numpy as np
import scipy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib as mpl

fal = 30 # fontsize axis labels
ftl = 20 # fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal

NS = 2 # system model degrees of freedom

def potential_energy(x, y, params_pe):
    """
    Potential energy function for the 2 DOF model
    x, y: N x N array of position values as meshgrid data
    params: 1 x 4 array of parameters
    """
    if np.size(x) == 1:
        nX = 1
        nY = 1
    else:
        nX = np.size(x,1)
        nY = np.size(x,0)
        x = np.ravel(x, order = 'C')
        y = np.ravel(y, order = 'C')
    
    # epsilon, Dx, alpha, lambd  = params_pe
    
    v1 = 0
    for i in range(5):
        v1 += params_pe[i]*(x**i) 
        
    v2 = params_pe[5]*((params_pe[6] - y)**2)

    vint = params_pe[7]/((y - x)**12)

    pe = np.reshape(v1 + v2 + vint, (nY, nX))
#     print(np.size(vy))
#     print(x,y,vyx)
    
    if np.size(x) == 1:
        pe = pe[0,0]

    return pe


def grad_hamiltonian(states, *params):

    # print(params[0][0:2]) 
    mu1 = params[0][0]
    mu2 = params[0][1]
    coeff_pe = params[0][2:]
    q1, q2, p1, p2 = states

    dvint_dq = (12*coeff_pe[-1])/((q2 - q1)**13)
    
    q1Dot = p1/mu1

    q2Dot = p2/mu2

    v1 = 0
    for i in range(5):
        v1 += i*coeff_pe[i]*(q1**(i-1))
    p1Dot = - ( v1 + dvint_dq )

    p2Dot = - ( - 2*coeff_pe[-3]*(coeff_pe[-2] - q2) \
                - dvint_dq )
    
    return np.array([q1Dot, q2Dot, p1Dot, p2Dot])


def jacobian(x, y, params):
    """
    Return the Jacobian of the Hamiltonian vector field for the solute-solvent model with Lennard-Jones term
    """

    mu1 = params[0]
    mu2 = params[1]
    coeff_pe = params[2:]

    jacobian_mat = np.zeros((2*NS,2*NS))

    d2vint = (156*coeff_pe[-1])/((y - x)**14)
    d2v1 = 0
    for i in range(5):
        d2v1 += i*(i - 1)*coeff_pe[i]*(x**(i - 2))
        
    deriv_wrt_r1_f3 = - d2v1 - d2vint
    
    deriv_wrt_r2_f3 = d2vint
    
    deriv_wrt_r1_f4 = d2vint
    
    deriv_wrt_r2_f4 = - 2*coeff_pe[-3] - d2vint
    
    
    jacobian_mat[2,0] = deriv_wrt_r1_f3

    jacobian_mat[2,1] = deriv_wrt_r2_f3
        
    jacobian_mat[3,0] = deriv_wrt_r1_f4
    
    jacobian_mat[3,1] = deriv_wrt_r2_f4

    
    jacobian_mat[0,2] = 1/mu1
    jacobian_mat[1,3] = 1/mu2

    return jacobian_mat



def vector_field(params, t, states):

    # print(params[0][0:2]) 
    mu1 = params[0]
    mu2 = params[1]
    coeff_pe = params[2:]
    q1, q2, p1, p2 = states

    dvint_dq = (12*coeff_pe[-1])/((q2 - q1)**13)
    
    q1Dot = p1/mu1

    q2Dot = p2/mu2

    v1 = 0
    for i in range(5):
        v1 += i*coeff_pe[i]*(q1**(i-1))
    p1Dot = - ( v1 + dvint_dq )

    p2Dot = - ( - 2*coeff_pe[-3]*(coeff_pe[-2] - q2) \
                - dvint_dq )
    
    return np.array([q1Dot, q2Dot, p1Dot, p2Dot])


def kinetic_energy(p1, p2, params_mass):
    """
    Returns the kinetic energy of the trajectory with same length 
    as the vector
    """
    ke = (1/2)*(p1**2/params_mass[0] + p2**2/params_mass[1])

    return ke
