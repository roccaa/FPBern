
% Time

% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32;

pi = ((( (-2/1) * x5 + ( - x7 + ( - x9 + ( - x10 + ( - x11 + ( - x12 + ( - x13 + ( - x14 + ( - x19 + ( - x24 + ( - x27 + ( - x29 + ( - x31 + ( - x32)))))))))))))) * x4) * x1 + ((( (2/1) * x5 + (x7 + (x9 + (x11 + (x12 + (x13 + (x14 + (x19 + (x24 + (x27 + (x29 + (x31 + (x32))))))))))))) * x4 + (x15 + (x16 + (x17 + (x18 + (x19 + (x24 + (x27 + (x29 + ( - x30)))))))))) * x2 + (((x5 + (x6 + (x7 + (x9 + (x12 + (x13 + (x14 + (x19 + (x24 + (x27 + (x29 + (x31 + (x32))))))))))))) * x4 + (x20 + (x21 + (x22 + (x23 + (x24 + (x27 + ( - x28)))))))) * x3 + (( - x5 + ( (-2/1) * x7 + ( - x9 + ( - x13 + ( - x14 + ( - x19 + ( - x24 + ( - x27 + ( - x29 + ( - x31 + ( - x32))))))))))) * x4^2)))) * x1 + ((( (-2/1) * x5 + ( - x15 + ( - x16 + ( - x17 + ( - x18 + ( - x19 + ( - x24 + ( - x27 + ( - x29 + ( - x31 + ( - x32))))))))))) * x2 + ((( - x5 + ( - x6 + ( - x7 + ( - x25 + ( - x26 + ( - x27 + ( - x29 + ( - x31 + ( - x32))))))))) * x4 + ( (2/1) * x5 + ( (2/1) * x6 + (x16 + (x17 + (x18 + (x19 + (x20 + (x21 + (x22 + (x23 + ( (2/1) * x24 + ( (2/1) * x27 + ( (2/1) * x29 + ( (2/1) * x31 + ( (2/1) * x32)))))))))))))))) * x3 + ((x5 + (x7 + (x17 + (x18 + (x19 + (x24 + (x27 + (x29 + (x31 + (x32)))))))))) * x4))) * x2 + ((( (-2/1) * x6 + ( - x21 + ( - x22 + ( - x23 + ( - x24 + ( - x27 + ( - x29 + ( - x31 + ( - x32))))))))) * x3 + ((x6 + (x7 + (x22 + (x23 + (x24 + (x27 + (x29 + (x31 + (x32))))))))) * x4)) * x3 + (( - x7 + ( - x32)) * x4)));
params = [x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32];
vars = [x1 x2 x3 x4];
X = [4,6.36,4,6.36,4,6.36,4,6.36];
tic
[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');
time = toc
disp(' ');
disp('The Result of Bernstein Computation of kepler1 with real number approx is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
res = max(abs(bmin),abs(bmax));
disp(' ');