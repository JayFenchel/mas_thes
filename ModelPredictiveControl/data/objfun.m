function [ f_val ] = objfun( U )
%OBJFUN Summary of this function goes here
%   Detailed explanation goes here
load data_matrix.mat


f_val = 0
for i=1:5
    x = U(30*(i)+1:30*(i+1));
    u = U(180+i);
    f_val = f_val + x'*Q_total*x + u'*R_total*u

end

