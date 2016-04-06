# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 22:12:31 2016

@author: LENOVO
"""

from scipy import exp,arange,array
from pylab import plot,show,title,xlabel,legend
from scipy.integrate import odeint

# Staem Methane Reaction CH4  +  2H2O  -->  CO2  +  4H2
# Constants

Ta      = 1000.0
H       =-40000.0
Cpm     = 40.0
Cmo     = 0.5
Fao     = 5.0
s       = 5.67*(10^(-8))
e       = 2/3
#Initial Conditions
f_int = [0.0,450.0,1.0]
w = arange(0,10,1)
#  Differential Equations
def pfr(f,w):
    x,T,y = f
    
    k  = 0.5*exp((41800.0/8.314)*(1.0/450.0-1.0/T))
    Cm = Cmo*(1-x)*450.0/T/(1.0+e*x)*y
    Cc = Cmo*0.5*x*(450.0/T)/(1.0+e*x)*y
    Kc = 165740.0*exp(H/8.314*(1.0/450.0-1.0/T))
    rA = -k*(Cm**2-Cc/Kc)
    
    dxdw = -rA/Fao
    dTdw = (0.8*(1*(Ta-T)+rA*H))/(Cpm*Fao) 
    dydw = -0.015*(1-0.5*x)*(T/450.0)/(2.0*y)
    
    return array([dxdw,dTdw,dydw], dtype=float) 
# Solve differential euations
fs = odeint(pfr,f_int,w)
# Plot the solution 
xs = fs[:,0]
Ts = fs[:,1]
ys = fs[:,2]
# Pull Concentration out for plotting
Cm = Cmo*(1-xs)*450.0/Ts/(1.0+e*xs)*ys
Cc = Cmo*0.5*xs*(450.0/Ts)/(1.0+e*xs)*ys
plot(w,xs,'-',w,Ts/1000.0,'--',w,ys,'-.',w,Cm,'o',w,Cc,'x')
title('Reaction Conversion, Temperature and concentration')
xlabel('L, Reactor length, m')
legend(['x', 'T,deg K','y','Cm','Cc'], 'best')
show()