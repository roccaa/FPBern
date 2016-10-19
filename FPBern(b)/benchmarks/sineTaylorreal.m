
% Time

% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14;

pi = (((( (-1/720) * x2 + ( (-1/5040) * x3 + ( (-1/5040) * x4 + ( (-1/5040) * x7 + ( (-1/5040) * x8 + ( (-1/5040) * x11 + ( (-1/5040) * x12 + ( (-1/5040) * x13 + ( (-1/5040) * x14))))))))) * x1^2 + ( (1/24) * x2 + ( (1/120) * x3 + ( (1/120) * x4 + ( (1/120) * x7 + ( (1/120) * x8 + ( (1/120) * x9 + ( (1/120) * x10 + ( (1/120) * x14))))))))) * x1^2 + ( (-1/2) * x2 + ( (-1/6) * x3 + ( (-1/6) * x4 + ( (-1/6) * x5 + ( (-1/6) * x6 + ( (-1/6) * x10 + ( (-1/6) * x14)))))))) * x1^2 + (x2 + (x6 + (x10 + (x14))))) * x1;
params = [x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14];
vars = [x1];
X = [-1.557079632679,1.557079632679,];
tic
[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');
time = toc
disp(' ');
disp('The Result of Bernstein Computation of sineTaylor is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
res = max(abs(bmin),abs(bmax));
disp(' ');