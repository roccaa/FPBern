
% Time
tic
% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27;

pi = (( (-2/1) * x7 + ( - x20 + ( - x21 + ( - x22 + ( - x23 + ( - x24 + ( - x25 + ( - x26 + ( - x27))))))))) * x1 + (( (2/1) * x7 + (x21 + (x22 + (x23 + (x24 + (x25 + (x26 + (x27)))))))) * x2 + ((x7 + (x8 + (x22 + (x23 + (x24 + (x25 + (x26 + (x27)))))))) * x3 + (( - x7 + ( - x9 + ( - x23 + ( - x24 + ( - x25 + ( - x26 + ( - x27))))))) * x4 + ((x7 + (x10 + (x24 + (x25 + (x26 + (x27)))))) * x5 + ((x7 + (x11 + (x25 + (x26 + (x27))))) * x6)))))) * x1 + ((( - x7 + ( - x8 + ( - x16 + ( - x17 + ( - x19 + ( - x27)))))) * x3 + ((x7 + (x10 + (x13 + (x15 + (x17 + (x19 + (x27))))))) * x5)) * x2 + (((x8 + (x11 + (x14 + (x15 + (x17 + (x19 + (x27))))))) * x6) * x3 + ((( - x10 + ( - x11 + ( - x18 + ( - x19 + ( - x27))))) * x6) * x5)));
params = [x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27];
vars = [x1 x2 x3 x4 x5 x6];
X = [4,6.36,4,6.36,4,6.36,4,6.36,4,6.36,4,6.36];

[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');

disp(' ');
disp('The Result of Bernstein Computation of kepler0 is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
toc
disp(' ');