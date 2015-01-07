#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 10:48:23 2014

@author: jayf
"""
from __future__ import division
import math
import threading
import numpy as np
import scipy.linalg

import rospy

from std_msgs.msg import *
from sensor_msgs.msg import Joy
from ottocar_msgs.msg import perception_laneStatePoly
from scipy.interpolate import interp1d

class Node(object):
    
    def __init__(self):
        super(Node, self).__init__()

        rospy.init_node('lane_controller')
        
        self.lock = threading.Lock()

        self.init_params()

        self.pub_angle_cmd = rospy.Publisher('/actuators/angle_cmd', Int8, tcp_nodelay=True,queue_size=1)
        self.pub_vel_cmd = rospy.Publisher('/functions/velocity_cmd', Float32, tcp_nodelay=True,queue_size=1)
        
        rospy.Subscriber('/joy', Joy, self.joy_callback)
        rospy.Subscriber('/functions/ottocar_perception/laneStatePoly', perception_laneStatePoly, self.callback_lanestatepoly)
        rospy.Subscriber('/functions/ottocar_perception/laneState', Float64MultiArray, self.callback_lanestate)
        rospy.Subscriber('/GPIO_button2', Bool, self.callback_button)
        
        self.run_status = False
        self.disable = False
        
        self.r = rospy.Rate(self.rate)
        
        
    def init_params(self):
        self.lanestate = None
        self.pos = None
        self.angle = None
        self.rate = 100
        
        
    def spin(self):
        rospy.loginfo("lane_controller_lqr started")
        while(not rospy.is_shutdown() and True):  # and self.run_status
            # print("Test klappt.")

            if (self.pos!=None) and (self.angle!=None):

                # Startwerte TO DO aktuelle Werte verwenden
                u = np.array([0., 0.])
                # TO DO Winkelgeschwindigkeit, Geschwindigkeit, Lenkeinschlag fehlen
                x = np.array([0., self.pos, self.angle, 0., 0., 1, 0.])

                L_R, L_G = self.modell_matrizen(x ,u)

                # print str(L_R) + "\n\n" + str(L_G)

                u[0] = np.dot(-L_R, np.matrix(([x[1]], [x[2]], [x[3]], [x[6]])))
                u[1] = np.dot(L_G, np.matrix(([abs(x[3])], [1 - x[5]])))

                # Stellgroessenbeschraenkung
                if u[0] > 0.35:
                    u[0] = 0.35
                    # print "positiv"
                elif u[0] < -0.35:
                    u[0] = float(-0.35)

                if u[1] > 2:
                    u[1] = 2
                elif u[1] < 0:
                    u[1] = 0

                print "\n\n" + str(u)


            self.r.sleep()

    # Modellbeschreibung
    def modell_matrizen(self, x, u):
        # Modellparameter
        lv = 0.1285
        lh = 0.1285
        m = 2.
        r = 0.0335
        g = 9.80665
        l = lv + lh
        rol = 1.2041
        T = 0.01
        gamma = 0.5
        epsilon = 8.78
        cwA = 1.15*0.025944
        fR = 4.*0.45
        Jz = 1/12 * m * (0.38 ** 2 + 0.185 ** 2)
        Cv = 1.
        Ch = 0.9

        Weg = 0

        # Systembeschreibung fuer den Einschlag der Raeder...
        A_R1 = [0, x[5] * math.cos(x[2] + x[4]), 0, 0]

        A_R2 = [0, 0, 1, 0]

        A_R3 = [0, 0, -((Ch*lh ** 2)/(x[5]*((lh ** 2*x[3] ** 2)/x[5] ** 2 + 1)) +
                        (Cv*lv ** 2*math.cos(x[6]))/(x[5]*((lv ** 2*x[3] ** 2)/x[5] ** 2 + 1)))/Jz,
                (Cv*lv*math.cos(x[6]) - Cv*lv*math.sin(x[6])*(x[6] - math.atan((lv*x[3])/x[5])) +
                 (gamma*lv*u[1]*math.cos(x[6]))/(epsilon*r))/Jz]

        A_R4 = [0, 0, 0, -1/T]

        A_R = np.matrix((A_R1, A_R2, A_R3, A_R4))

        B_R = np.matrix(([0], [0], [0], [1/T]))  # Spaltenvektor

        # ...fuer die Geschwindigkeit
        A_G1 = [A_R3[2], ((Ch*lh ** 2*x[3])/(x[5] ** 2*((lh ** 2*x[3] ** 2)/x[5] ** 2 + 1)) +
                          (Cv*lv ** 2*x[3]*math.cos(x[6]))/(x[5] ** 2*((lv ** 2*x[3] ** 2)/x[5] ** 2 + 1)))/Jz]

        if (Weg >= np.spacing(1e10) and x[6] >= np.spacing(1e13)) or (Weg >= np.spacing(1e8) and x[6] >= np.spacing(1e15)):
            A_G2 = [((Ch*lh*math.sin(x[4]))/(x[5]*((lh**2*x[3]**2)/x[5]**2 + 1)) -
                     (Cv*lv*math.sin(x[4] - x[6]))/(x[5]*((lv**2*x[3]**2)/x[5]**2 + 1)) +
                     (Weg*m*x[5]**2*math.cos(x[4])*math.acos(1 - (Weg**2*math.tan(x[6])**2)/(2*lh**2)) *
                      ((lh*lv*math.cos(x[4] - x[6] + (lv*x[3])/x[5]))/(l*x[5]) -
                       (lh*lv*math.cos(x[4] - (lh*x[3])/x[5]))/(l*x[5])))/g)/m,
                    -(fR + cwA*rol*x[5] - (Cv*lv*x[3]*math.sin(x[4] - x[6]))/(x[5]**2*((lv**2*x[3]**2)/x[5]**2 + 1)) +
                      (Ch*lh*x[3]*math.sin(x[4]))/(x[5]**2*((lh**2*x[3]**2)/x[5]**2 + 1)) +
                      (Weg*m*x[5]**2*math.cos(x[4])*math.acos(1 - (Weg**2*math.tan(x[6])**2)/(2*lh**2)) *
                       ((lh*lv*x[3]*math.cos(x[4] - x[6] + (lv*x[3])/x[5]))/(l*x[5]**2) -
                        (lh*lv*x[3]*math.cos(x[4] - (lh*x[3])/x[5]))/(l*x[5]**2)))/g -
                      (2*Weg*m*x[5]*math.cos(x[4])*math.acos(1 - (Weg**2*math.tan(x[6])**2)/(2*lh**2)) *
                       ((lh*math.sin(x[4] - x[6] + (lv*x[3])/x[5]))/l + (lv*math.sin(x[4] - (lh*x[3])/x[5]))/l))/g)/m]
        else:
            A_G2 = [((Ch*lh*math.sin(x[4]))/(x[5]*((lh**2*x[3]**2)/x[5]**2 + 1)) -
                     (Cv*lv*math.sin(x[4] - x[6]))/(x[5]*((lv**2*x[3]**2)/x[5]**2 + 1)))/m,
                    -(fR + cwA*rol*x[5] - (Cv*lv*x[3]*math.sin(x[4] - x[6]))/(x[5]**2*((lv**2*x[3]**2)/x[5]**2 + 1)) +
                      (Ch*lh*x[3]*math.sin(x[4]))/(x[5]**2*((lh**2*x[3]**2)/x[5]**2 + 1)))/m]
        A_G = np.matrix((A_G1, A_G2))
        B_G = np.matrix(([(gamma*lv*math.sin(x[6]))/(Jz*epsilon*r)],
                         [((gamma*math.cos(x[4] - x[6]))/(epsilon*r) - (math.cos(x[4])*(gamma - 1))/(epsilon*r))/m]))

        # LQR-Regler: Gewichtungsmatrizen Q und R
        Q_R = np.diag([1., 0.1, 0.001, 1e-6])
        R_R = np.matrix(1e-6)   # Die Typen von SciPy sind Idioten..., nur deswegen.

        L_R = self.lqr(A_R, B_R, Q_R, R_R)

        Q_G = np.diag([1, 0.0001])
        R_G = np.matrix(1e-5)

        L_G = self.lqr(A_G, B_G, Q_G, R_G)

        return L_R, L_G

    # Berechnung der lqr-optimalen Zustandsrückführung
    def lqr(self, a, b, q, rr):
        # aus: https://github.com/markwmuller/controlpy
        xx = scipy.linalg.solve_continuous_are(a, b, q, rr)
        ll = np.dot(scipy.linalg.inv(rr), (np.dot(b.T, xx)))
        return ll

    def lqr_discrete_time(self, a, b, q, rr):
        # aus: https://github.com/markwmuller/controlpy
        xx = scipy.linalg.solve_discrete_are(a, b, q, rr)
        ll = np.dot(scipy.linalg.inv(np.dot(np.dot(b.T, xx), b)+rr), (np.dot(np.dot(b.T, xx), a)))
        return ll

    def callback_lanestatepoly(self,msg):
        if(not self.run_status):
            return

        self.lock.acquire()
        
        self.time_alt = self.time
        self.time = rospy.get_time() ####msg.header.ZEIT
        self.poly_confidence_last = self.poly_confidence
        self.poly_confidence = msg.confidence.data

        if self.poly_confidence >= 0.5:
            self.poly_last_time = self.time

        xs = np.array(msg.targetPolyXArray_VRF.data)
        ys = np.array(msg.targetPolyYArray_VRF.data)
        
        if xs.shape[0]>=2:
            s=np.zeros_like(xs)
            for i in range(1, size(xs)):
                s[i] = s[i-1] + np.sqrt(np.square(xs[i]-xs[i-1]) + np.square(ys[i]-ys[i-1]))
            
            if xs.shape[0]<3:
                print "linear (<3)"
                fx=interp1d(s, xs, kind='linear')
                fy=interp1d(s, ys, kind='linear')
            else:
                try:
                    fx=interp1d(s, xs, kind='linear')
                    fy=interp1d(s, ys, kind='linear')
                    print "cubic possible"
                except:
                    print "linear"
                    fx=interp1d(s, xs, kind='linear')
                    fy=interp1d(s, ys, kind='linear')
                    
        else:
            rospy.loginfo("Keine verwendbare neue Polyline")
                        

        self.lock.release()

    def callback_lanestate(self,msg):
        # example msg
        # layout:
        #   dim: []
        #   data_offset: 0
        # data: [946685062.8521258, 120.0, -0.0, 0.5]

        self.lanestate = msg.data
        self.pos = self.lanestate[1]            # in cm
        self.angle = self.lanestate[2]          # in grad
        
        
    def callback_button(self,msg):
        self.run_status = not self.run_status
        self.start_zeit = rospy.get_time()
    
    
    def joy_callback(self, joy):
        if(not self.run_status):
            return
        self.joy = joy
        if self.joy.buttons[5]==1:  # RB pushed
            self.disable = True
        else:
            self.disable = False
        
        
if __name__ == '__main__':
    process = Node()
    process.spin()