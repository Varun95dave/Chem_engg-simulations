# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 21:49:00 2016

@author: vdave
"""
from math import factorial
from scipy import exp
import scipy
from scipy.optimize import fmin
from scipy.integrate import trapz
from pylab import plot,show
C=scipy.array([0,2,3,2,1.5,1.2,0])
t=scipy.array([0,1,2,3,4,5,6])
I=trapz(C,t)
C=C/I
V=2.0
Q=4.0
tao=V/Q
def nonideal(n):
    taoi=tao/n
    e=(t**(n-1))/factorial(n-1)/taoi**(n-1)*exp(-t/taoi)
    I=trapz(e,t)
    e=e/I
    e=(e-C)**2
    e=sum(e)
    return e
    
n=fmin(nonideal,3)
print n
taoi=tao/n
e=(t**(n-1))/factorial(n-1)/taoi**(n-1)*exp(-t/taoi)
I=trapz(e,t)
e=e/I
plot(t,C)
plot(t,e)
show()