# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 22:12:31 2016

@author: vdave
"""

from scipy import exp,arange,array
from pylab import plot,show,title,xlabel,legend
from scipy.integrate import odeint

# Reaction 2M  -->  C
# Constants
Q       = 10 #Volumetric flow rate
Ta      = 2000.0  # Surrounding temp (Flame temp)
H       =-40000.0 # Heat of reaction (exothermic)
Fmo     = 5.0 # Molar flow rate
Cmo     = Fmo/Q # Pure feed to the pfr
s       = 5.67*(10^(-8)) # sigma: stefan-boltzma's constant
e       = 1/2 # epsilon
k0      = 0.5 #rate constant at tref
Kc0     = 165740.0 #eqm constant at tref

# Constants for Heat capacity of product [1]
A       = 7
B       = 10
C       = -4
D       = 5
E       = 8

#Initial Conditions
f_int = [0.0,450.0,1.0]
w = arange(0,15,1)
#  Differential Equations
def pfr(f,w):
    x,T,y = f
    
    # Variation of rate constant (K) with temperature
    k  = k0*exp((41800.0/8.314)*(1.0/450.0-1.0/T))
    
    # Variation of reactant (Cm) with temperature
    Cm = Cmo*(1-x)*450.0/T/(1.0+e*x)*y
    
    # Variation of product (Cc) with temperature
    Cc = Cmo*0.5*x*(450.0/T)/(1.0+e*x)*y
    
    # Variation of equilibrium constant (Kc) with temperature
    Kc = Kc0*exp(H/8.314*(1.0/450.0-1.0/T))
    
    # Rate law
    rA = -k*(Cm**2-Cc/Kc)
    dxdw = -rA/Fmo
    
    #Heat Capcity of reactant varying with temp
    Cpm=A+B*T+C*T*T+D*T*T*T+E/(T*T)
    
    #Temperature gradient across pfr length
    a=T*T*T*T
    b=Ta*Ta*Ta*Ta
    dTdw = (0.8*(100*(Ta-T)+rA*H)+s*(b-a))/(Cpm*Fmo)  #Accounting for radiative heat Loss
    
    #Pressure gradient across pfr length
    dydw = -0.015*(1-0.5*x)*(T/450.0)/(2.0*y)
    
    return array([dxdw,dTdw,dydw], dtype=float) 
# Solve differential euations
fs = odeint(pfr,f_int,w)

# Plot the solution 
xs = -fs[:,0]
Ts = fs[:,1]
ys = fs[:,2]

# Plot
plot(w,xs,'-',w,Ts/100.0,'--',w,ys,'-.')
title('Reaction Conversion, Temperature and pressure')
xlabel('x, Length of pfr, m')
legend(['x', 'T,deg K','y'], 'best')
show()
