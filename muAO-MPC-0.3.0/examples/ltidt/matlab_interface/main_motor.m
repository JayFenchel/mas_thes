%% Basic example
% This file creates the system module as a mat file.
% It then calls the (extenal) Python code generator.
% Finally compile and use the generated MPC controller.
%% create the interface
create_sys_motor_mat();
!python main_motor.py
%% compile it and make it available
cd cmpc/matlab
make
cd ../..
addpath ./cmpc/matlab
%% Usage example
ctl = mpc_ctl;  % create an instance of the class
ctl.conf.in_iter = 10;  % configure the algorithm
x = [0.1; -0.5];  % current state
% forming a QP is only necessary if a different QP solver is used
ctl.form_qp(x)  % form the QP for the current state
qpx = ctl.qpx;  % qpx contains the QP in standard form
u = quadprog(qpx.HoL, qpx.gxoL, [], [], [], [], qpx.u_lb, qpx.u_ub);
% the QP is automatically formed when using its own algorithm (ALM+FGM) 
ctl.solve_problem(x);  % solve the MPC problem for current state
ctl.u_opt  % display the computed input sequence
norm(u - ctl.u_opt)  % ctl.u_opt is the MPC approximation of u
