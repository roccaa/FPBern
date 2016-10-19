
% Time

% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4;

pi = (x1*x1-x1)*(x4)+x1*x1*x3+(2*x1*x1-x1)*x2;
params = [x2 x3 x4];
vars = [x1];
X = [0,1];

tic
[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');
time = toc
disp(' ');
disp('The Result of Bernstein Computation for simple_test is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
res = max(abs(bmin),abs(bmax));

disp(' ');