# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 15:17:56 2016

@author: vdave
"""

import scipy
from scipy.integrate import odeint

import matplotlib.pyplot as plt
#for the absorber
Gin=100#Input Moles of Gas
GCin=0.5*Gin#Input Moles of CO2
GMin=0.5*Gin#Input Moles of CH4
GWin=0#Input Moles of Water Vapor
Lin=100#Input Moles of Liquid
LCin=0#Input Moles of CO2 in liquid
LMin=0#Input Moles of CH4 in Liquid
LWin=1*Lin#Input Moles of Liquid Water

def eqn(F,z):
    #Data
    KgC=0.1 #in kmol/hr.m3.kPa(Mass Transfer Coefficient of CO2)
    KgM=0.1 #in kmol/hr.m3.kPa(Mass Transfer Coefficient of CH4)
    KgW=0.1 ##in kmol/hr.m3.kPa(Mass Transfer Coefficient of Water)
    KlC=0.1
    KlM=0.1
    KlW=0.1
    HC1=0.034#Henry's Constant for CO2(in mol/lit.atm)
    HM1=0.009#Henry's Constant for CH4(in mol/lit.atm)
    HC=0.034*1.013
    HM=0.009*1.013
    HW=1.013
    Psat=3.170 # psat of water @ 298K in kPa
    P=1 # Operating pressure in atm
    Cc=0#Concentration of CO2 in Liquid
    Cm=0#Concentration of CH4 in liquid
    Cw=1#Concentration of water
    h=10#Height of Column
    A=10#Area of Cross section
    n=10#No of Divisions 
    gc1=0
    gm1=0
    gw1=0
    dz=float(h/(n-1))#Length of one division
    dgcdz=-1*KgC*A*(F[0]/(F[0]+F[1]+F[2])*P-Cc/HC)#Solving all the ODEs
    dgmdz=-1*KgM*A*(F[1]/(F[0]+F[1]+F[2])*P-Cm/HM)
    dgwdz=-1*KgW*A*(F[2]/(F[0]+F[1]+F[2])*P-Psat)
    dlcdz=-1*KlC*A*(F[3]/(F[3]+F[4]+F[5])*HC*P-(Cc))
    dlmdz=-1*KlM*A*(F[4]/(F[3]+F[4]+F[5])*HM*P-(Cm))
    dlwdz=-1*KlW*A*(F[5]/(F[3]+F[4]+F[5])*HW*P-(Cw))
    Cc=F[3]/(F[5]/55.55)#IN Mol/lit
    Cm=F[4]/(F[5]/55.55)
    Cw=F[5]/(F[5]/55.55)
    
    return(dgcdz,dgmdz,dgwdz,dlcdz,dlmdz,dlwdz)#Returns the values of the ODE Solved
    
z=scipy.arange(0,10,1)#Array to Store the X-axis values

F0=[50,50,0,0,0,100]#Initial Values

value=odeint(eqn,F0,z)#Calling odeint to solve the set of Equations

FgC=value[:,0]#Stores the values of CO2 in gas at various intervals
FgM=value[:,1]#Stores the values of CH4 in gas at various intervals
FgW=value[:,2]#Stores the values of water in gas at various intervals
FlC=value[:,3]#Stores the values of CO2 in liquid at various intervals
FlM=value[:,4]#Stores the values of CH4 in liquid at various intervals
FlW=value[:,5]#Stores the values of water in liquid at various intervals

e=(FlW[0]-FlW[9])/FlW[0]#Error in Liquid water calculations

print ("CO2 in gas: ",FgC) #Printing the values 
print ("CH4 in gas: ",FgM)
print ("Water in gas",FgW)
print ("CO2 in liquid",FlC)
print ("CH4 in liquid",FlM)
print ("Liquid water",FlW)
print("Error in Liquid water calculations : ",e)

plt.plot(z,FgC)#Plotting the values of the moles
plt.plot(z,FgM)
plt.plot(z,FgW)
plt.plot(z,FlC)
plt.plot(z,FlM)
plt.plot(z,FlW)
plt.xlabel('Height in m')
plt.ylabel('Flow Rate in mol/hr')
plt.show()
'''References:
   1.Perry's Chemical Engineering Handbook-1999(pg no 2114 & 2115)
   2.Henry's Constants: Henry-3.0 pdf,R.Sander:Henry's Law Constant
   3.Antoine Co-efficients:
       http://www.eng.auburn.edu/~drmills/mans486/Diffusion%Tube/Antoine_coefficient_table.PDF
   '''
