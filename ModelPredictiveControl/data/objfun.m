function [ f_val ] = objfun( U )
%OBJFUN Summary of this function goes here
%   Detailed explanation goes here
load data_matrix.mat
ref = zeros(30, 1);
ref(19) = 200;

% f_val = U(31:60)'*Q_total*U(31:60)...
%     + U(61:90)'*Q_total*U(61:90)...
%     + U(91:120)'*Q_total*U(91:120)...
%     + U(121:150)'*Q_total*U(121:150)...
%     + U(151:180)'*Q_total*U(151:180)...
%     + U(181)'*R_total*U(181)...
%     + U(182)'*R_total*U(182)...
%     + U(183)'*R_total*U(183)...
%     + U(184)'*R_total*U(184)...
%     + U(185)'*R_total*U(185);

f_val = 0;
for i=1:5
    x = U(30*(i)+1:30*(i+1));
    u = U(180+i);
    f_val = f_val + (x-ref)'*Q_total*(x-ref) + u'*R_total*u;
end

end

