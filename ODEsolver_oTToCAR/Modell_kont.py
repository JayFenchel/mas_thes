#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 16:00:24 2014

@author: jayf
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.integrate import odeint

from time import sleep, time

x0, t0 = np.array([0, 0, 0, 0, 0, 1, 0]), 0
t1 = 1
dt = .01
n_sample = (t1-t0)/dt + 1
dim_x = x0.shape[0]


class OdeModel:

    def __init__(self):
        # Physikalische Eigenschaften
        # Gemessen
        m = 2.          # [kg] Gewicht                                          *nicht ueberprueft*
        r = 0.0335      # [m] Radradius                                         *nicht ueberprueft*
        g = 9.80665     # [m/s^2] Erdbeschleunigung
        l = 0.2570      # [m] Gesamtlaenge des Radstandes                       *nicht ueberprueft*
        rho_l = 1.2041  # [kg/m^3] Luftdichte (Meeresspiegel, 20 Grad Celsius)

        # Aus Ueberlegung geschaetzt
        m_vertlg = 0.5          # Verhaeltnis lv/l (lv - Abstand Vorderachse zum Fahrzeugschwerpunkt)
        T = 0.2                 # [s] zeitliche Verzoegerung des Inputs
        J = 1./12. * m *(np.square(0.38) + np.square(0.185)) # [kg*m^2] Traegheitsmoment eines Qauders mit den Kantenlaengen a=0.38m, b=0.185m
        c_wA = 1.15*0.025944    # [-]*[m^2] Stroemungswiderstandskoeffizient mal Flaeche

        # Aus Messdaten gschaetzt
        gamma = 0.5061
        xi = 84.5
        f_r = 10.6577
        c_wert = 24.4686
        c_alpha_vertlg = 0.3844

        xi_v = xi
        xi_h = xi
        c_h_alpha = c_wert

        # Zusammengesetzte Parameter für Zustandsgleichungen
        self.lv=l*m_vertlg;    # [m] Laenge vordere Haelfte
        self.lh=l-self.lv;            # [m] Laenge hintere Haelfte
        
        c_v_alpha=c_h_alpha*c_alpha_vertlg;
        
        self.p1=c_h_alpha*self.lh/J;
        self.p2=c_v_alpha*self.lv/J;
        self.p3=gamma*self.lv/J/r/xi_v;
        self.p4=c_v_alpha/m;
        self.p5=gamma/m/r/xi_v;
        self.p6=c_h_alpha/m;
        self.p7=(1.-gamma)/r/xi_h/m;
        self.p8=f_r/m;
        self.p9=(c_h_alpha*np.square(m/2.)+c_v_alpha*np.square(m/2.))/m/c_v_alpha/c_h_alpha;
        self.p10=c_wA*rho_l/2./m;
        self.p11=T;


    def f(self, x, t, u):

        u=np.array([0.1, 30])
    ## einfacheres Testsystem
        # dx1 = x[4]*np.cos(x[2])
        # dx2 = x[4]*np.sin(x[2])
        # dx3 = x[3]
        # dx4 = 0
        # dx6 = 0
        #
        # return [dx1, dx2, dx3, dx4, dx6]

    ## Meine Modellgleichungen aus Matlab
        dx0 = x[5]*np.cos(x[2]+x[4])
        dx1 = x[5]*np.sin(x[2]+x[4])
        dx2 = x[3]
        dx3 = (-self.p1*np.arctan(self.lh*x[3]/x[5]) + self.p2*x[6]*np.cos(x[6])
            - self.p2*np.arctan(self.lv*x[3]/x[5])*np.cos(x[6]) + self.p3*u[1]*x[6])
        dx4 = (-x[3] + self.p4*x[6]*np.cos(x[6]-x[4])/x[5] - self.p4*np.arctan(self.lv*x[3]/x[5])*np.cos(x[6]-x[4])/x[5]
            + self.p5*x[6]*u[1]/x[5] - self.p5*x[4]*u[1]/x[5] + self.p6*np.arctan(self.lh*x[3]/x[5])/x[5]
            - self.p7*x[4]*u[1]/x[5] + np.square(self.p9*x[3])*x[4]*x[5])
        dx5 = (self.p4*x[6]*x[4] - np.square(self.p4*x[6]) - self.p4*np.arctan(self.lv*x[3]/x[5])*x[6]
            + self.p4*np.arctan(self.lv*x[3]/x[5])*x[4] + self.p5*u[1]*np.cos(x[6]-x[4])
            + self.p6*np.arctan(self.lh*x[3]/x[5])*x[4] + self.p7*u[1] - self.p8*x[5]
            - np.square(self.p9*x[5])*np.square(x[3]) - self.p10*np.square(x[5]))
        # dx5 = 0
        dx6 = -x[6]/self.p11 + u[0]/self.p11

    #    dx0 = x[5]*np.cos(x[2]+x[4]);
    #    dx1 = x[5]*np.sin(x[2]+x[4]);
    #    dx2 = x[3];
    #    dx3 = 1 / Jz * (Fvq * lv * cos(x[6]) + Fvl * lv * sin(x[6]) - Fhq * lh);
    #    dx4 = (1 / (m * x[5])) * (Fvq * cos(x[6]-x[4]) + Fvl * sin (x[6]-x[4]) + Fhq * cos(x[4]) + FK * sin(x[4]) - Fhl * sin(x[4])) - x[3];
    #    dx5 = (1/m) * (Fvl * cos(x[6] - x[4]) + Fhq * sin(x[4]) + Fhl * cos(x[4]) - Fvq * sin(x[6]-x[4]) - FR - FL - FK * cos(x[4]));
    #    dx6 = -x[6]/T + u0/T;
    #    return [dx0, dx1, dx2, dx3, dx4, dx5, dx6]
        return [dx0, dx1, dx2, dx3, dx4, dx5, dx6]

myDGL = OdeModel()
t = np.linspace(0, 1, 21)
y0=x0
time0 = time()   

for i in range(20):
    soln = odeint(myDGL.f, y0, t, args = (0,))

print(time()-time0)
plt.figure(1)
plt.subplot(3,1,1)
plt.title('kontinuierliches System')
plt.plot(soln.T[0],soln.T[1])
plt.subplot(3,1,2)
plt.plot(t,soln.T[2])
plt.subplot(3,1,3)
plt.plot(t,soln.T[5])
plt.show()
# input()
