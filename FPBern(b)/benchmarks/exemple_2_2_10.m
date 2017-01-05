% Time

% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16;
pi = (( 22 * x3 + ( 22 * x5 + ( 11 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11 + ( 7 * x12 + ( 8 * x13 + ( 9 * x14 + ( 10 * x15 + ( 11 * x16))))))))))))) * x1 + (( 22 * x3 + ( 22 * x4 + ( 44 * x5 + ( 22 * x6 + ( 4 * x7 + ( 6 * x8 + ( 8 * x9 + ( 10 * x10 + ( 12 * x11 + ( 14 * x12 + ( 16 * x13 + ( 18 * x14 + ( 20 * x15 + ( 22 * x16)))))))))))))) * x2)) * x1 + (( 22 * x4 + ( 22 * x5 + ( 11 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11 + ( 7 * x12 + ( 8 * x13 + ( 9 * x14 + ( 10 * x15 + ( 11 * x16))))))))))))) * x2^2);
params = [x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16];
vars = [x1 x2];
X = [-1,1,-1,1];
tic
[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');
time = toc
disp(' ');
disp('The Result of Bernstein Computation of exemple-2-2-10 with real number approx is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
res = max(abs(bmin),abs(bmax));
disp(' ');