%% System module example as MATLAB mat file
% Based on the motor system described in the basic tutorial.
% To run the complete example, run the main_motor.m script
clear
%% auxiliary constanst
T = 0.1;
K = 0.2;
dt = T / 10;
N = 10;
%% continuous time system
Ac = [0, 1; 0, -1/T];
Bc = [0; K/T];
% input constraints
u_lb = -100;
u_ub = 100;
%% weighting matrices
Q = [1, 0; 0, 1];
R = 1;
P = 'auto'; 
%% system module
save sys_motor.mat *
