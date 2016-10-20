
% Time

% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18;

pi = (( (4/1) * x4 + ( (2/1) * x7 + (x8 + (x11 + ( (2/1) * x13 + (x14 + ( (2/1) * x15 + ( (2/1) * x18)))))))) * x1 + (( (-2/1) * x4 + ( (-2/1) * x5 + ( (-4/1) * x7 + ( (-2/1) * x8 + ( (-2/1) * x11 + ( (-2/1) * x15 + ( (-2/1) * x18))))))) * x2 + (( (-2/1) * x4 + ( (-4/1) * x6 + ( (-2/1) * x12 + ( (-4/1) * x13 + ( (-2/1) * x14 + ( (-2/1) * x15 + ( (-2/1) * x18))))))) * x3^2))) * x1 + ((( (4/1) * x5 + ( (2/1) * x7 + (x8 + ( (2/1) * x9 + (x10 + ( (2/1) * x11 + ( (2/1) * x15 + ( (2/1) * x18)))))))) * x2 + ( (-2/1) * x5 + ( (-4/1) * x9 + ( (-2/1) * x10 + ( (-2/1) * x11 + ( (-2/1) * x15 + ( (-2/1) * x18))))))) * x2 + (((( (4/1) * x6 + ( (2/1) * x12 + ( (2/1) * x13 + (x14 + (x15 + (x18)))))) * x3^2 + ( (2/1) * x6 + ( (2/1) * x16 + (x17 + (x18))))) * x3 + ( (-2/1) * x6 + ( (-4/1) * x16 + ( (-2/1) * x17 + ( (-2/1) * x18))))) * x3 + ( (2/1) * x9 + (x10 + (x11 + (x15 + ( (2/1) * x16 + (x17 + ( (2/1) * x18)))))))));
params = [x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18];
vars = [x1 x2 x3];
X = [-10,10,-10,10,-10,10];
tic
[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');
time = toc
disp(' ');
disp('The Result of Bernstein Computation of schwefel is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
res = max(abs(bmin),abs(bmax));
disp(' ');
