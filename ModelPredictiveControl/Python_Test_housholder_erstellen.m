% Solves a LLS problem min || A x - b || with Householder
% Sebastian Sager, 2013

J1 = [1 -2 0 .1; 1 -4 0 0; 1 -6 1 0; 1 -8 0 1];
F1 = [4.5; 2.5; 2; 0.99];

A = J1;
b = F1;

[neta,np] = size(A);

Q = eye(neta);
R = A;
% Householder Spalten
for i=1:np
    y = R(i:neta,i);
    alpha = norm(y);
    v = (y - alpha*eye(neta-i+1,1));
    v = v / norm(v);

    H = eye(neta,neta);
    H(i:neta,i:neta) = eye(neta-i+1) - 2*v*v';
    Q = H*Q;
    R = H*R;
end
Q
R
% Aufloesung
c = Q*b;
Rtilde = R(1:np,1:np);
ctilde = c(1:np);
x = Rtilde \ ctilde;
residual = 0;	
for i=np+1:neta
    residual = residual + 0.5*c(i)^2;
end
