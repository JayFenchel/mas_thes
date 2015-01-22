#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 16:00:24 2014

@author: jayf
"""
import numpy as np


class DiscreteOdeModel:

    def __init__(self):
        # Physikalische Eigenschaften
        # Gemessen
        m = 2.          # [kg] Gewicht                      *nicht ueberprueft*
        r = 0.0335      # [m] Radradius                     *nicht ueberprueft*
        g = 9.80665     # [m/s^2] Erdbeschleunigung
        l = 0.2570      # [m] Gesamtlaenge des Radstandes   *nicht ueberprueft*
        rho_l = 1.2041  # [kg/m^3] Luftdichte (Meeresspiegel, 20 Grad Celsius)
        # Aus Ueberlegung geschaetzt
        m_distribution = 0.5  # Verhaeltnis lv/l (lv - Abstand Vorderachse zum Fahrzeugschwerpunkt)
        T = 0.2               # [s] zeitliche Verzoegerung des Inputs
        J = (1./12.*m*(np.square(0.38) + np.square(0.185)))  # [kg*m^2] Traegheitsmoment eines Qauders mit den Kantenlaengen a=0.38m, b=0.185m
        c_wA = 1.15*0.025944  # [-]*[m^2] Stroemungswiderstandskoeffizient mal Flaeche

        # Aus Messdaten gschaetzt
        gamma = 0.5061
        xi = 84.5
        f_r = 10.6577
        c_value = 24.4686
        c_alpha_distribution = 0.3844

        xi_v = xi
        xi_h = xi
        c_h_alpha = c_value

        c_v_alpha=c_h_alpha*c_alpha_distribution

        # Zusammengesetzte Parameter für Zustandsgleichungen
        self.lv=l*m_distribution  # [m] Laenge vordere Haelfte
        self.lh=l-self.lv         # [m] Laenge hintere Haelfte

        self.p1 = c_h_alpha*self.lh/J
        self.p2 = c_v_alpha*self.lv/J
        self.p3 = gamma*self.lv/J/r/xi_v
        self.p4 = c_v_alpha/m
        self.p5 = gamma/m/r/xi_v
        self.p6 = c_h_alpha/m
        self.p7 = (1.-gamma)/r/xi_h/m
        self.p8 = f_r/m
        self.p9 = ((c_h_alpha*np.square(m/2.) + c_v_alpha*np.square(m/2.))
                / m/c_v_alpha/c_h_alpha)
        self.p10 = c_wA*rho_l/2./m
        self.p11 = T

    def f(self, x, delta_t, u):

        ## Meine Modellgleichungen aus Matlab
        dx0 = (x[5]*np.cos(x[2]+x[4]))*delta_t + x[0]
        dx1 = (x[5]*np.sin(x[2]+x[4]))*delta_t + x[1]
        dx2 = (x[3])*delta_t + x[2]
        dx3 = (-self.p1*np.arctan(self.lh*x[3]/x[5]) + self.p2*x[6]*np.cos(x[6])
            - self.p2*np.arctan(self.lv*x[3]/x[5])*np.cos(x[6]) + self.p3*u[1]*x[6])*delta_t + x[3]
        dx4 = (-x[3] + self.p4*x[6]*np.cos(x[6]-x[4])/x[5] - self.p4*np.arctan(self.lv*x[3]/x[5])*np.cos(x[6]-x[4])/x[5]
            + self.p5*x[6]*u[1]/x[5] - self.p5*x[4]*u[1]/x[5] + self.p6*np.arctan(self.lh*x[3]/x[5])/x[5]
            - self.p7*x[4]*u[1]/x[5] + np.square(self.p9*x[3])*x[4]*x[5])*delta_t + x[4]
        dx5 = (self.p4*x[6]*x[4] - np.square(self.p4*x[6]) - self.p4*np.arctan(self.lv*x[3]/x[5])*x[6]
            + self.p4*np.arctan(self.lv*x[3]/x[5])*x[4] + self.p5*u[1]*np.cos(x[6]-x[4])
            + self.p6*np.arctan(self.lh*x[3]/x[5])*x[4] + self.p7*u[1] - self.p8*x[5]
            - np.square(self.p9*x[5])*np.square(x[3]) - self.p10*np.square(x[5]))*delta_t + x[5]
        # dx5 = 0
        dx6 = (-x[6]/self.p11 + u[0]/self.p11)*delta_t + x[6]

        return [dx0, dx1, dx2, dx3, dx4, dx5, dx6]

    def f_t_first(self, t, x, u):

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

        return [dx0, dx1, dx2, dx3, dx4, dx5, dx6]