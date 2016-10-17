
% Time
tic
% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13;

pi = (( (-2/1) * x4 + ( - x7 + ( - x8 + ( - x11 + ( - x12 + ( - x13)))))) * x2 + ( - x4 + ( - x12 + ( - x13)))) * x1 + ((( (-2/1) * x4 + ( (-2/1) * x5 + ( (-2/1) * x9 + ( (-2/1) * x10 + ( (-2/1) * x11 + ( (-2/1) * x12 + ( (-2/1) * x13))))))) * x3) * x2 + (( - x5 + ( - x13)) * x3));
params = [x4 x5 x6 x7 x8 x9 x10 x11 x12 x13];
vars = [x1 x2 x3];
X = [-15,15,-15,15,-15,15];

[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');

disp(' ');
disp('The Result of Bernstein Computation of rigidBody1 is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
toc
disp(' ');