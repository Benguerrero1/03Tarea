# -*- coding: utf-8 -*-
"""
Created on Wed Oct 07 14:31:16 2015

@author: Administrador
"""

'''
Script que realiza RK3 para integrar la ecuacion de van der Pol
y grafica el resultado.
'''

import numpy as np
import matplotlib.pyplot as plt

mu=1.370

def f(y,v):
    return (v,-y-mu*(y**2-1)*v)

def get_k1(y_n,v_n,h,f):
    f_n = f(y_n,v_n)
    return h*f_n[0],h*f_n[1]

def get_k2(y_n,v_n,h,f):
    k1 = get_k1(y_n,v_n,h,f)
    f_n = f(y_n+k1[0]/2.,v_n+k1[1]/2.)
    return k1,(h*f_n[0],h*f_n[1])

def get_k3(y_n,v_n,h,f):
    k1,k2 = get_k2(y_n,v_n,h,f)
    f_n = f(y_n-k1[0]+2*k2[0],v_n-k1[1]+2*k2[1])
    return k1,k2,(h*f_n[0],h*f_n[1])

def RK3 (y_n,v_n,h,f):
    '''
    Recibe los valores en el paso n-esimo de "y" y "v"
    y retorna los valores en el paso siguiente.
    '''
    k1,k2,k3 = get_k3(y_n,v_n,h,f)
    y_n1 = y_n + 1./6. * (k1[0]+4*k2[0]+k3[0])
    v_n1 = v_n + 1./6. * (k1[1]+4*k2[1]+k3[1])
    return y_n1,v_n1
    