#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 10:48:23 2014

@author: jayf
"""
import math
import numpy as np

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
        rospy.Subscriber('/GPIO_button2', Bool, self.callback_button)
        
        self.run_status = False
        self.disable = False
        
        self.r = rospy.Rate(self.rate)
        
        
    def init_params(self):
        self.rate = 100
        
        
    def spin(self):
        rospy.loginfo("lane_controller_lqr started")
        while(not rospy.is_shutdown() and self.run_status):
            
            
    def optimal_control
        
        #SOLL OPTIMAL CONTROL PROBLEM LÖSEN
    
    
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