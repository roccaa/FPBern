
% Time
tic
% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16;

pi = (((( (-5/32) * x2 + ( (-1/32) * x5 + ( (-5/128) * x12 + ( (-5/128) * x13 + ( (-5/128) * x14 + ( (-5/128) * x15 + ( (-5/128) * x16))))))) * x1 + ( (3/16) * x2 + ( (1/16) * x4 + ( (1/16) * x8 + ( (1/16) * x9 + ( (1/16) * x10 + ( (1/16) * x11 + ( (1/16) * x16)))))))) * x1 + ( (-1/4) * x2 + ( (-1/8) * x3 + ( (-1/8) * x5 + ( (-1/8) * x6 + ( (-1/8) * x7 + ( (-1/8) * x11 + ( (-1/8) * x16)))))))) * x1 + (x2 + ( (1/2) * x3 + ( (3/2) * x4 + ( (3/2) * x7 + ( (3/2) * x11 + ( (3/2) * x16))))))) * x1 + (x4 + (x7 + (x11 + (x16))));
params = [x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16];
vars = [x1];
X = [0,1];

[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');

disp(' ');
disp('The Result of Bernstein Computation of sqroot is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
toc
disp(' ');