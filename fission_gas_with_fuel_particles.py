    # -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 15:19:22 2017
@author: Parth
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 17:10:21 2017
@author: parth
"""

import numpy as np
import math
import time
import os
from matplotlib import pyplot as plt

number = int(1E6)
dt = 1E-5

#% Initialisations
T = np.zeros(number)
Db = np.zeros(number)
A = np.zeros(number)
v = np.zeros(number) #bubble rising velocity m/s
t = np.zeros(number)
y = np.zeros(number)
rhol = np.zeros(number)
rhog = np.zeros(number)
Pb = np.zeros(number)
N = np.zeros(number) #chack
Ma = np.zeros(number)
N_aerosol = np.zeros(number)
Mg = np.zeros(number)
par = np.zeros(number)
RF = np.zeros(number)
rho_aero = np.zeros(number)
Vb = np.zeros(number)
diff = np.zeros(number)
sedimentation = np.zeros(number)
inertial_impaction = np.zeros(number)
number_dens_removed = np.zeros(number)
#%% Constants
#mass_of_vap_sod = 1E-3 # mass of vaporised sodium in kg

Tb = 1154.1 # Boiling point of sodium in kelvine
Ts = 931 # sodium pool temperature
cp = 1284.4 #kJ/kg
#lg = 4197E3
k = 1.38E-23 #Boltzman constant J/K
g = 9.8 # gravitational acceleration
P0 = 4E8 # Atmosperic pressure 1 bar = 1E5 pascals
H = 6.234 # height of the pool
h = 5 #check
k = 1.3807E-23 #boltzman constant in J/k
#dp = 0.1E-6 #particle diameters
Na = 6.022E23 #avogrado number
#lembda = 0.693/(33120)
#lembda = 0.693/1E10
lembda = 1e-2
#Va = 1/6 * math.pi * dp**3
molar_mass_gas = 135E-3 #(xe-135)
molar_mass_aerosol = 235E-3 #(fuel U-235)
c = 0
#rhop = 4510
rhop =1e4 #input('enter the value of rhop:  ')
#% Initial Values


#Pb[0] = 2E5 #check
T[0]= 661 +273+1000 #check
rhol[0]=949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2 #density of sodium
Pb[0] = P0+ rhol[0]*g*(H)
v[0] =0.01 #check value
N[0] = 2.1136E24
rhog[0] = (N[0]*molar_mass_gas/Na)

Db[0]= (6*N[0]*k*T[0]/(3.14*Pb[0]))**(1/3.0)


mug = 2.28E-5 #gas viscosity (xenon)



rho_aero[0] = rhop

print('sodium density \t{:0.3e}'.format(rhol[0]))
print('bubble diameter \t{:0.3e}'.format(Db[0]*1E2))
aerosol_size  =[1e-8,1e-7,1e-6,1e-5]
#aerosol_size = np.linspace(0.01E-6,10e-6,int(10))
#aerosol_size = [1E-6,1E-5]
DF = np.zeros([len(aerosol_size),number])
j = 0
print(Pb[0])
Vb[0] = math.pi * 0.1667 * Db[0]**3
Ma[0] = rho_aero[0]*Vb[0] #initial mass of aerosols in bubble
#%%  Program loop
for dp in aerosol_size:
    print(dp)
    for i in range(number-3):

        t[i] = i*dt

        y[i] = (i)*dt*v[i]

        Ma[i] = Ma[0]*math.exp(-y[i]*c)

        rhol[i] = 949-0.223*(T[i]-273.15)-1.7E-5*(T[i]-273.15)**2 #density of sodium

        Pb[i] =P0 + rhol[i]*g*(H-y[i]) # pressure equation

        if i > 0 :
            Db[i] = (6*N[i]*k*T[i]/(3.14*Pb[i]))**(1/3.0) # diameter equation

        Vb[i] = math.pi*0.1667*Db[i]**3

        rho_aero[i] = Ma[i]/(0.1667*math.pi*Db[i]**3)

        rhog[i] =  (N[i]*molar_mass_gas/Na)

#        number_dens_removed[i] = (- (rho_aero[i]-rho_aero[i-1])*Na /   molar_mass_aerosol)

        if i>0:
            number_dens_removed[i] = (Ma[i]-Ma[i-1])*Na / 0.235

        N[i+1] = N[i] + dt*(-lembda*N[i]) + number_dens_removed[i] *  dt

        T[i+1] = T[i] + dt* (0-h*3.14*Db[i]**2 * (T[i]-Ts))/(rhog[i]*Vb[i]*cp) # temperature equation


        v[i+1] = v[i]+ dt*(6/(rhol[i]*3.14*Db[i]**3))*(  (3.14*Db[i]**3 *g*(rhol[i]-rhog[i]))/6 - 12*3.14*Db[i]*v[i]) #velocity of bubble

        lembd = k*T[i]/(float(math.sqrt(2))*math.pi*Pb[i]*Db[i]**2) #mean free path of a gas molecule

        cunn = 1 + (2*lembd/dp) * (1.257 + 0.4 * math.exp(-0.55*dp/lembd)) #cunningham slip correction

        theta = k*T[i]*cunn/(3*3.14*mug*dp)

        tau = rho_aero[i]*dp**2 *cunn /(18*mug)

        diff[i] = 1.8*(8*theta/(v[i]*Db[i]**3))**(1/2.0) #diffusion coefficient

        sedimentation[i] = 1.5*g*tau/(Db[i]*v[i]) #sedimentation coefficient

        inertial_impaction[i] = 18*v[i]*tau/Db[i]**2 # inertial impaction coefficient

        c = diff[i] + sedimentation[i] + inertial_impaction[i]

        dy = dt*v[i]
        par[i] = dy*(-c*Ma[i])
        Ma[i+1] = Ma[i] + par[i]

        RF[i] = math.exp(-y[i]*c)


    #    N[i+1] = rhog[i+1]*Na/235
    #    if i >0:
    #        Ma[i] =  Ma[0] * math.exp(-c*y[i])

        if y[i]>H:
            print('bubble reached the pool surface')
            break
        elif N[i]<0:
            print('bubble gas decayed out')
            break

#        plt.savefig('Release')

    DF[j,0:i] = RF[0:i]

    j = j+1

for j in range(len(aerosol_size)):
    plt.plot(t[0:i-1],DF[j,0:i-1],label='{:0.1e}'.format(aerosol_size[j]))
plt.legend()
plt.xlabel('time(second)')
plt.ylabel('fraction in bubble')
plt.savefig('fraction of mass in bubble')
'''
plt.figure()
plt.plot(aerosol_size,DF[:,i-1])
plt.semilogx()
plt.semilogy()
'''
'''
plt.plot(aerosol_size*1e6,DF[:j])
plt.semilogx()
plt.semilogy()

plt.figure()
plt.plot(t[0:i],Db[0:i],label='diameter',color='blue')
plt.xlabel('time(seconds)',color='blue',size=14)
plt.ylabel('diameter (meter)',color='blue',size=14)
plt.legend()
plt.savefig('diameter')
'''


'''
#plt.figure()
plt.plot([i for i in range(0,j)],(DF[0:j]),label='release fraction')

plt.xlabel('Time',color='blue',size=14)
plt.ylabel('fraction in bubble',color='blue',size=14)
plt.legend()
plt.savefig('Release')
'''
'''
plt.figure()
plt.plot(t[:i],T[:i],label = 'Temperature',color='blue')
plt.xlabel('time',color='blue',size=14)
plt.ylabel('Temp (K)',color='blue',size=14)
plt.legend()
plt.savefig('temperature')
plt.figure()
plt.plot(t[:i],Pb[:i],label='pressure',color='blue')
plt.xlabel('time(seconds)',color='blue',size=14)
plt.ylabel('Pressure (Pa)',color='blue',size=14)
plt.legend()
plt.savefig('pressure')
plt.figure()
plt.plot(t[:i],v[:i],label='velocity',color='blue')
plt.xlabel('time(second)',color='blue',size=14)
plt.ylabel('velocity(m/s)',color='blue',size=14)

plt.legend()
plt.savefig('velocity')

plt.figure()
plt.plot(t[:i],N[:i],label='number density',color='blue')
plt.xlabel('time(second)',color='blue',size=14)
plt.ylabel('number density',color='blue',size=14)

plt.legend()
plt.savefig('number density')

'''

'''
#convert meter to cm
Db =(Db*1E2)
y = y*1E2-250
#convert array to list
Db = list(Db)
y = list(y)
import turtle

def draw_circle(turtle, color, size, x, y):
    turtle.penup()
    turtle.color(color)

    turtle.goto(x,y)
    turtle.pendown()

    turtle.circle(size)


turtle.shape("circle")
#turtle.resizemode("user")
turtle.shapesize(1E-2,1E-2)
#tommy.speed(7)
scale = 100
# Draw three circles:
for i in range(int(len(y)/scale)):
    draw_circle(turtle, "green", Db[scale*i], 25, y[scale*i])



'''



