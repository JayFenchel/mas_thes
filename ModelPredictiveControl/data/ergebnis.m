close all
clear all

load myproblem_bounds.mat
load data_matrix.mat

problem = optimproblem
problem.options.MaxIter = 20

steps = 25
x = zeros(185, steps);
for i=1:steps
    step = i
    [x(:, i)] = fmincon(problem);
    problem.x0(1:30) = Asys*problem.x0(1:30) + Bsys*problem.x0(181);
    problem.x0(31:150) = x(61:180, i);
    problem.x0(151:180) = x(151:180, i);
    problem.x0(181:184) = x(182:185, i);
    problem.x0(185) = x(185, i);
    problem.lb(1:30) = problem.x0(1:30);
    problem.ub(1:30) = problem.x0(1:30);
end


E_matlab = x(19, :)

% Python Algorithmus (05.06.2015, 25 steps)
E_python = [4.00000000e+02   3.99959534e+02   3.98559660e+02   3.85976302e+02...
    3.56261174e+02   3.14042535e+02   2.67653506e+02   2.24785668e+02...
    1.91094408e+02   1.69524307e+02   1.60062806e+02   1.60594853e+02...
    1.68223676e+02   1.79737206e+02   1.91986728e+02   2.02545217e+02...
    2.09996536e+02   2.13916804e+02   2.14653990e+02   2.12942539e+02...
    2.09749667e+02   2.06118288e+02   2.02867540e+02   2.00471870e+02...
    1.99099548e+02]

figure(1)
    hold on
    grid on
    plot(E_matlab)
    plot(E_python)