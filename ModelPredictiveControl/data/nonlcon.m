function [ c, ceq ] = nonlcon( U )
%OBJFUN Summary of this function goes here
%   Detailed explanation goes here
load socp_matrices.mat

c(1) = norm(V21*U) + M21*U + c_max;
c(2) = norm(V22*U) + M22*U + c_max;
c(3) = norm(V23*U) + M23*U + c_max;
c(4) = norm(V24*U) + M24*U + c_max;
c(5) = norm(V25*U) + M25*U + c_max;

ceq = 0;

end