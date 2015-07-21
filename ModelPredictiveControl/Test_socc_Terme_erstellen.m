clc
format rat

T = 3;
n = 2;
m = 1;

z = [1 0 0 -3 -1 1 7 1 6]';

A = [1 7; 3 -4];
b = [11; 5];
c = [3; 7];
d = -1;
AE = [.5 7; 3 -4];
bE = [10.5; 5];
cE = [3; 6.5];
dE = -.5;

zi = z(8:9)
if ((A*zi + b)'*(A*zi + b) - (norm(A*zi+b))^2) > 0.00000000001
    'Warnung'
end
%% P_of_zk
nab1_f = -2*((cE'*zi + dE)*cE - AE'*(AE*zi + bE))';

% P_ref = np.array([[0., 58., 128., 0., 0., 0., 0., 0., 0.],
%                   [0., 0., 0., 0., 4., 212., 0., 0., 0.],
%                   [0., 0., 0., 0., 0., 0., 0., -292., 661./2.]])

%% form_d
minus_f = (cE'*zi + dE)^2 - (AE*zi + bE)'*(AE*zi + bE);
d_i= 1/minus_f

% d_ref = np.array([[-1./145.],
%                   [-1./284.],
%                   [-4./5371.]])

%% Term in Phi    
minus_f = (c'*zi + d)^2 - (A*zi + b)'*(A*zi + b);
nab2_f = -2*(c*c' - A'*A);

Teil_Term = 1/minus_f*nab2_f;

% term_for_socc_ref = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
%                               [0, -2/145, 52/145, 0, 0, 0, 0, 0, 0],
%                               [0, 52/145, -32/145, 0, 0, 0, 0, 0, 0],
%                               [0, 0, 0, 0, 0, 0, 0, 0, 0],
%                               [0, 0, 0, 0, -1/142, 13/71, 0, 0, 0],
%                               [0, 0, 0, 0, 13/71, -8/71, 0, 0, 0],
%                               [0, 0, 0, 0, 0, 0, 0, 0, 0],
%                               [0, 0, 0, 0, 0, 0, 0, -2/5371, 224/5371],
%                               [0, 0, 0, 0, 0, 0, 0, 224/5371, -182/5371]])