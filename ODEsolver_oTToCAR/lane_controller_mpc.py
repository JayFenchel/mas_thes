#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 10:48:23 2014

@author: jayf
"""
import math
import threading
import numpy as np

import rospy

from std_msgs.msg import *
from sensor_msgs.msg import Joy
from ottocar_msgs.msg import perception_laneStatePoly
from ottocar_msgs.msg import StampedFloat32
from ottocar_controller import DiscreteOdeModel
from scipy.interpolate import interp1d
from scipy.optimize import fmin

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
        rospy.Subscriber('/GPIO_button2', Bool, self.callback_button)
        rospy.Subscriber('/functions/yaw_rate', StampedFloat32, self.callback_yaw_rate)
        rospy.Subscriber('/functions/velocity', StampedFloat32, self.callback_velocity)
        
        self.run_status = False
        self.disable = False
        
        self.r = rospy.Rate(self.rate)
        
        
    def init_params(self):
        self.yaw_rate = None
        self.velocity = None
        self.fx = None
        self.fy = None
        self.u_last = 0  # erster Startwert for Optimierung

        self.gradeaus = - 30  # geschaetzt am 08.01.2015

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
        rospy.loginfo("lane_controller_mpc started")
        while(not rospy.is_shutdown() and self.run_status):
            print("Test klappt.")

            if ((self.fx != None) and (self.fy != None) and
                    (self.yaw_rate != None) and (self.velocity != None)):

                x0 = np.array([0., self.fx, self.fy, self.yaw_rate, 0., self.velocity, 0.])

                u = self.optimal_control(x0)

                u[0] = min(max(u[0], -127), 127)
                self.u_last = u  # letzter Wert als Startwert for Optimierung
                angle_cmd = min(max(u[0] + self.gradeaus, -127), 127)

                self.pub_angle_cmd.publish(angle_cmd)

            self.r.sleep()

    def cost_function_diskret(self, u1, xzero):
        # TODO Sollwerte aus Polyline übergeben
        discrete_model = DiscreteOdeModel()  # TODO Nicht jedes Mal neu initiieren
        dt = 0.1  # TODO global definieren
        us = [u1, 30]
        for i in range(1,5,1):
            solution = np.array(discrete_model.f(xzero, dt/5, us))
            xzero = solution
        return np.square(solution.T[2]-0).sum() + .01*np.square(u1)

    def optimal_control(self, x0):
        x0 = x0  # np.array([0, 0, 0.2, 0, 0, 1, 0])
        u0 = self.u_last  # [0.1, 30]
        #SOLL OPTIMAL CONTROL PROBLEM LÖSEN
        uopt = fmin(self.cost_function_diskret, u0[0], args=(x0,))

        return uopt
    
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
                self.fx=interp1d(s, xs, kind='linear')
                self.fy=interp1d(s, ys, kind='linear')
            else:
                try:
                    self.fx=interp1d(s, xs, kind='linear')
                    self.fy=interp1d(s, ys, kind='linear')
                    print "cubic possible"
                except:
                    print "linear"
                    self.fx=interp1d(s, xs, kind='linear')
                    self.fy=interp1d(s, ys, kind='linear')
                    
        else:
            rospy.loginfo("Keine verwendbare neue Polyline")

        self.lock.release()

    def callback_yaw_rate(self, msg):
        self.yaw_rate = msg.float32.data            # in rad

    def callback_velocity(self, msg):
        self.velocity = msg.float32.data            # in m/s
        
        
    def callback_button(self,msg):
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
    process.spin()