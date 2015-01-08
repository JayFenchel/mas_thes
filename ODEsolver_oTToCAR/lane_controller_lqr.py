#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 10:48:23 2014

@author: jayf
"""
from __future__ import division
import math
import numpy as np
import scipy.linalg

import rospy

from std_msgs.msg import Bool
from std_msgs.msg import Int8
from std_msgs.msg import Float32
from std_msgs.msg import Float64MultiArray
from sensor_msgs.msg import Joy
from ottocar_msgs.msg import perception_laneStatePoly
from ottocar_msgs.msg import StampedFloat32

class Node(object):
    
    def __init__(self):
        super(Node, self).__init__()

        rospy.init_node('lane_controller')

        self.init_params()

        self.pub_angle_cmd = rospy.Publisher('/actuators/angle_cmd', Int8, tcp_nodelay=True,queue_size=1)
        self.pub_vel_cmd = rospy.Publisher('/functions/velocity_cmd', Float32, tcp_nodelay=True,queue_size=1)
        
        rospy.Subscriber('/joy', Joy, self.joy_callback)
        rospy.Subscriber('/functions/ottocar_perception/laneState', Float64MultiArray, self.callback_lanestate)
        rospy.Subscriber('/GPIO_button2', Bool, self.callback_button)
        rospy.Subscriber('/functions/yaw_rate', StampedFloat32, self.callback_yaw_rate)
        rospy.Subscriber('/functions/velocity', StampedFloat32, self.callback_velocity)

        self.run_status = False
        self.disable = False
        
        self.r = rospy.Rate(self.rate)

        self.spin_restarter()

    def init_params(self):
        self.lanestate = None
        self.pos = None
        self.angle = None
        self.lane_confidence = None
        self.lane_confidence_last = 0
        self.yaw_rate = None
        self.velocity = None
        self.angle_cmd_last = 0


        self.gradeaus = - 30  # geschaetzt am 08.01.2015
        self.glaettung = 0.8  # 1 - nur neuen Wert verwenden
        self.decrease = 150  # Dämpfung der Lenkstellgröße aus dem LQR
        self.rate = 100

    def spin_restarter(self):
        print "waiting for button"
        while(not rospy.is_shutdown()):
            if(self.run_status):
                print("Button gedrückt - LaneController läuft")
                self.spin()
                print("Button gedrückt - LaneController gestoppt")
            self.r.sleep()

    def spin(self):
        rospy.loginfo("lane_controller_lqr started")
        # Startwerte für Stellgrößen
        u = np.array([0., 1.])
        while(not rospy.is_shutdown() and self.run_status):  # and self.run_status
            # print("Test klappt.")

            if (self.pos != None) and (self.angle != None) and\
                    (self.yaw_rate != None) and (self.velocity != None):

                # TO DO Lenkeinschlag fehlt, überhaupt nötig?
                x = np.array([0., self.pos, self.angle, self.yaw_rate, 0., self.velocity, 0.])

                L_R, L_G = self.modell_matrizen(x, u)

                u[0] = np.dot(-L_R, np.matrix(([20 - x[1]], [x[2]], [x[3]], [x[6]])))
                u[0] /= -self.decrease  # "-", um falsches Vorzeichen zu korrigieren
                u[1] = np.dot(L_G, np.matrix(([abs(x[3])], [1 - x[5]])))
                u[1] = 1

                # Stellgroessenbeschraenkung

                u[0] = min(max(u[0], -127), 127)
                angle_cmd = min(max(u[0] + self.gradeaus, -127), 127)

                print "\n\n" + str(angle_cmd)

                angle_cmd = (self.glaettung*self.lane_confidence*angle_cmd
                                 + (1-self.glaettung)*self.lane_confidence_last*self.angle_cmd_last)
                if abs(self.lane_confidence) + abs(self.lane_confidence_last) != 0:
                    angle_cmd /= (self.glaettung*self.lane_confidence
                                      + (1-self.glaettung)*self.lane_confidence_last)

                # velocity_cmd oder speed_cmd erstmal noch mit joystick
                # publishen. Dazu 'roslaunch ottocar_launch NUC...' launchen.
                self.pub_angle_cmd.publish(angle_cmd)
                self.angle_cmd_last = angle_cmd


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

    def callback_lanestate(self, msg):
        # example msg
        # layout:
        #   dim: []
        #   data_offset: 0
        # data: [946685062.8521258, 120.0, -0.0, 0.5]

        self.lanestate = msg.data
        self.pos = self.lanestate[1]                # in cm
        self.angle = self.lanestate[2]*2*np.pi/360  # in rad
        self.lane_confidence_last = self.lane_confidence
        self.lane_confidence = self.lanestate[3]

    def callback_yaw_rate(self, msg):
        self.yaw_rate = msg.float32.data            # in rad

    def callback_velocity(self, msg):
        self.velocity = msg.float32.data            # in m/s
        
        
    def callback_button(self, msg):
        self.run_status = not self.run_status
    
    
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