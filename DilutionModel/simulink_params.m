%%% sim params %%%

alpha = 0.06;
s0 = 0.70;
m = 0.5;

chi = [ 0.055; 0.05 ];
s = s0 - m*chi;
nu = min(max(1 - chi/alpha, zeros(2,1)), ones(2,1));

% s = [ 0.68; 0.72 ];
% nu = [ 0.3; 0.8 ];

f_10 = 0.1;
f_0 = [ f_10; 1 - f_10 ];

