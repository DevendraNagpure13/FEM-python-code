# -*- coding: utf-8 -*-
"""
Created on Thu May  4 00:55:06 2023

@author: devna
"""

import numpy as np
from scipy.integrate import quad 
import matplotlib.pyplot as plt
n= int(input("Enter the value of n: "))
K= np.zeros((n,n))
F= np.zeros(n)
e= 2.718
for i in range(n):
    for j in range(n):
       #bi= lambda x: x**i
       #fi= lambda x: (x**i)*np.sin(x)
       fi= lambda x: (x**i)*(1/(1+np.exp(-x)))                          
       kij= lambda x: (x**i)*(x**j)
       K[i,j],err= quad(kij, -2*np.pi, 2*np.pi)
       F[i], errf= quad(fi, -2*np.pi, 2*np.pi)
       
c= np.linalg.solve(K, F)

x= np.linspace(-2*np.pi, 2*np.pi, 100)
y= 1/(1+np.exp(-x))

def product(c,x):
    yh= np.zeros(100)
    for i in range(100):
        for j in range(n):
            yh[i]= yh[i]+ c[j]*(x[i]**j)
    return yh
yh= product(c,x)
plt.plot(x, y, 'ko')
plt.plot(x, yh, 'r*')
plt.show()
#%%







