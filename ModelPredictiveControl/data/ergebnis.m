close all
clear all

% load myproblem_bounds.mat
% load data_matrix.mat
% 
% problem = optimproblem;
% 
% % % Für -400 als Startwert
% % problem.x0 = -1*problem.x0
% % problem.ub(19) = -1*problem.ub(19)
% % problem.lb(19) = -1*problem.lb(19)
% problem.options.MaxIter = 30;
% 
% steps = 2
% x = zeros(185, steps);
% for i=1:steps
%     step = i
%     [x(:, i)] = fmincon(problem);
%     problem.x0(1:30) = Asys*x(1:30, i) + Bsys*x(181, i);
%     problem.x0(31:150) = x(61:180, i);
%     problem.x0(151:180) = x(151:180, i);
%     problem.x0(181:184) = x(182:185, i);
%     problem.x0(185) = x(185, i);
%     problem.lb(1:30) = problem.x0(1:30);
%     problem.ub(1:30) = problem.x0(1:30);
% end
% 
% 
% E_matlab = x(19, :)

% Python Algorithmus (05.06.2015, mittags, 25 steps)
E_python = [4.00000000e+02   3.99959534e+02   3.98559660e+02   3.85976302e+02...
    3.56261174e+02   3.14042535e+02   2.67653506e+02   2.24785668e+02...
    1.91094408e+02   1.69524307e+02   1.60062806e+02   1.60594853e+02...
    1.68223676e+02   1.79737206e+02   1.91986728e+02   2.02545217e+02...
    2.09996536e+02   2.13916804e+02   2.14653990e+02   2.12942539e+02...
    2.09749667e+02   2.06118288e+02   2.02867540e+02   2.00471870e+02...
    1.99099548e+02]

E_python2 = [  4.00000000e+02   3.99550153e+02   3.88544475e+02   3.60211002e+02...
    3.18806246e+02   2.72479529e+02   2.28719061e+02   1.93150185e+02...
    1.68972812e+02   1.57010865e+02   1.55883385e+02   1.62626757e+02...
    1.73802990e+02   1.86379686e+02   1.98056344e+02   2.07269645e+02...
    2.13213894e+02   2.15798971e+02   2.15502143e+02   2.13145436e+02...
    2.09672378e+02   2.05956619e+02   2.02701519e+02   2.00341091e+02...
    1.99014695e+02]

E_python3 = [  4.00000000e+02   3.99558373e+02   3.88744415e+02   3.60706300e+02...
    3.19635751e+02   2.73598199e+02   2.29558054e+02   1.92927299e+02...
    1.67144749e+02   1.53708039e+02   1.51776724e+02   1.58450435e+02...
    1.70122261e+02   1.83434800e+02   1.95979246e+02   2.06187093e+02...
    2.13145046e+02   2.16663391e+02   2.17448000e+02   2.16339278e+02...
    2.13890894e+02   2.10613542e+02   2.07093115e+02   2.03894650e+02...
    2.01417933e+02]
% v shifted, letzten werte neu initialisiert
E_python4 = [  4.00000000e+02   3.99558373e+02   3.88743339e+02   3.60674664e+02...
    3.19438268e+02   2.73184597e+02   2.29397603e+02   1.93703759e+02...
    1.69741147e+02   1.58250819e+02   1.57429281e+02   1.64195963e+02...
    1.75236421e+02   1.87607948e+02   1.99019577e+02   2.07942596e+02...
    2.13616483e+02   2.15985427e+02   2.15538883e+02   2.13093042e+02...
    2.09578760e+02   2.05855350e+02   2.02612856e+02   2.00272963e+02...
    1.98966789e+02]
E_python5 = [  4.00000000e+02   3.99789698e+02   3.94399298e+02   3.75160420e+02...
    3.40551608e+02   2.96925592e+02   2.52441027e+02   2.13873441e+02...
    1.85361345e+02   1.68296012e+02   1.61838821e+02   1.63706433e+02...
    1.70976118e+02   1.80756927e+02   1.90651349e+02   1.98993871e+02...
    2.04895571e+02   2.08147400e+02   2.09041982e+02   2.08168603e+02...
    2.06223473e+02   2.03861689e+02   2.01602236e+02   1.99785082e+02...
    1.98571151e+02]

% 21.07.15 - socc-Terme korrigiert und getestet
E_python6 = [  4.00000000e+02   3.99792399e+02   3.94774624e+02   3.83171141e+02...
    3.67231895e+02   3.48770710e+02   3.28900549e+02   3.08279311e+02...
    2.87320592e+02   2.66382975e+02   2.46064988e+02   2.27433451e+02...
    2.11757939e+02   1.99997589e+02   1.92499864e+02   1.88981890e+02...
    1.88689762e+02   1.90625554e+02   1.93764689e+02   1.97220321e+02...
    2.00338735e+02   2.02728954e+02   2.04241198e+02   2.04913548e+02...
    2.04905639e+02]

E_C = [
        400.000000
        399.792475
        394.776450
        383.175433
        367.238048
        348.778262
        328.909168
        308.288759
        287.330712
        266.393598
        246.075201
        227.440485
        211.757449
        199.985330
        192.473369
        188.941425
        188.638415
        190.568591
        193.708443
        197.170944
        200.301086
        202.705842
        204.233097
        204.918782
        204.920936
        204.455568
        203.744001
        202.976272
        202.291132
        201.770441
];

figure(1)
    hold on
    grid on
%     plot(E_matlab)
    plot(E_python6, 'r')
%     plot(E_python3, 'g')
    plot(E_python5, 'g')
    plot(E_C, 'b')