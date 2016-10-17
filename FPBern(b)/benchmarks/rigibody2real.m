
% Time
tic
% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18;

pi = (((( (-3/1) * x4 + ( - x5 + ( - x13 + ( - x14 + ( - x15 + ( - x16 + ( - x17 + ( - x18)))))))) * x3) * x2 + (( (4/1) * x4 + ( (2/1) * x5 + ( (2/1) * x7 + ( (2/1) * x8 + ( (2/1) * x9 + ( (2/1) * x12 + ( (2/1) * x16 + ( (2/1) * x17 + ( (2/1) * x18))))))))) * x3)) * x2) * x1 + (( - x4 + ( - x18)) * x2 + (( (12/1) * x5 + ( (6/1) * x10 + ( (6/1) * x11 + ( (3/1) * x12 + ( (3/1) * x16 + ( (6/1) * x17 + ( (6/1) * x18))))))) * x3^2));
params = [x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18];
vars = [x1 x2 x3];
X = [-15,15,-15,15,-15,15];

[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');

disp(' ');
disp('The Result of Bernstein Computation of rigidBody2 is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
toc
disp(' ');