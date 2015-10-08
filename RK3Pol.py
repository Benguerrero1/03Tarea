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
    
'''
Para condiciones iniciales v0=0 y0=0.1
'''
v0 = 0
y0 = 0.1

n = 628
h = 0.1
y_1 = np.zeros(n)
v_1 = np.zeros(n)

y_1[0] = y0
v_1[0] = v0

for i in range(1,n):
    (y_1[i],v_1[i]) = RK3(y_1[i-1],v_1[i-1],h,f)

plt.figure(1,figsize=(7,13))
plt.clf
plt.subplot(311)
plt.plot(y_1,v_1,color="r",label="condiciones iniciales: dy/ds=0, y=0.1")
plt.xlabel('$y$', fontsize=20)
plt.ylabel("$\\frac{dy}{ds}$",fontsize=20)
plt.legend(loc='lower right',prop={'size':10})
plt.title("Trayectoria oscilador de Van der Pol")

'''
Para condiciones iniciales v0=0 y0=4
'''
v0 = 0
y0 = 4

n = 628
h = 0.1
y_2 = np.zeros(n)
v_2 = np.zeros(n)

y_2[0] = y0
v_2[0] = v0

for i in range(1,n):
    (y_2[i],v_2[i]) = RK3(y_2[i-1],v_2[i-1],h,f)

plt.subplot(312)
plt.plot(y_2,v_2,color="b",label="condiciones iniciales: dy/ds=0, y=4")
plt.xlabel('$y$', fontsize=20)
plt.ylabel("$\\frac{dy}{ds}$",fontsize=20)
plt.legend(loc='lower right',prop={'size':10})

plt.subplot(313)
plt.plot(y_1,v_1,color="r",label="condiciones iniciales: dy/ds=0, y=0.1")
plt.plot(y_2,v_2,color="b",label="condiciones iniciales: dy/ds=0, y=4")
plt.xlabel('$y$', fontsize=20)
plt.ylabel("$\\frac{dy}{ds}$",fontsize=20)
plt.legend(loc='lower right',prop={'size':10})
plt.savefig("vdpol.png")

s_values=np.linspace(1,n*h,n)

plt.figure(2,figsize=(7,7))
plt.subplot(211)
plt.title("y(s) para condiciones iniciales: dy/ds=0, y=0.1")
plt.plot(s_values,y_1,color="r")
plt.ylabel("y(s)")

plt.subplot(212)
plt.title("y(s) para condiciones iniciales: dy/ds=0, y=4")
plt.plot(s_values,y_2,color="b")
plt.xlabel("s")
plt.ylabel("y(s)")
plt.savefig("vdpol2.png")
plt.show()