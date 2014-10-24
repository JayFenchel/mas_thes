%% Basic Simulink example
% This file creates the code for sample model and then opens simulink
% system, which allows to simulate it's dynamics.
%% create requirements file for controller
create_sys_motor_mat();
%% generate controller code
!python main_motor.py
%% open the Simulink model
open('simulink_example');
open_system('simulink_example/S-Function Builder')
