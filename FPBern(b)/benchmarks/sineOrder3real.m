
% Time
tic
% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7;

pi = ((( (-5120184131/10000000000) * x2 + ( (-1290061377/10000000000) * x4 + ( (-1290061377/10000000000) * x5 + ( (-1290061377/10000000000) * x6 + ( (-1290061377/10000000000) * x7))))) * x1 + ( (1/2) * x2 + ( (1/2) * x3 + ( (1/2) * x7)))) * x1 + ( (4774648293/5000000000) * x2 + ( (4774648293/5000000000) * x3 + ( (4774648293/5000000000) * x7)))) * x1
params = [x2 x3 x4 x5 x6 x7];
vars = [x1];
X = [-2,2];

[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');

disp(' ');
disp('The Result of Bernstein Computation of sineOrder3 is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
toc
disp(' ');

