% Time

% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11;
pi = ((((( 15 * x3 + ( 15 * x5 + ( 3 * x6 + ( 3 * x7 + ( 3 * x8 + ( 3 * x9 + ( 2 * x10 + ( 3 * x11)))))))) * x1 + (( 60 * x3 + ( 15 * x4 + ( 75 * x5 + ( 15 * x6 + ( 15 * x7 + ( 15 * x8 + ( 15 * x9 + ( 10 * x10 + ( 15 * x11))))))))) * x2)) * x1 + (( 90 * x3 + ( 60 * x4 + ( 150 * x5 + ( 30 * x6 + ( 30 * x7 + ( 30 * x8 + ( 30 * x9 + ( 20 * x10 + ( 30 * x11))))))))) * x2^2)) * x1 + (( 60 * x3 + ( 90 * x4 + ( 150 * x5 + ( 30 * x6 + ( 30 * x7 + ( 30 * x8 + ( 30 * x9 + ( 20 * x10 + ( 30 * x11))))))))) * x2^3)) * x1 + (( 15 * x3 + ( 60 * x4 + ( 75 * x5 + ( 15 * x6 + ( 15 * x7 + ( 15 * x8 + ( 15 * x9 + ( 10 * x10 + ( 15 * x11))))))))) * x2^4)) * x1 + (( 15 * x4 + ( 15 * x5 + ( 3 * x6 + ( 3 * x7 + ( 3 * x8 + ( 3 * x9 + ( 2 * x10 + ( 3 * x11)))))))) * x2^5);
params = [x3 x4 x5 x6 x7 x8 x9 x10 x11];
vars = [x1 x2];
X = [-1,1,-1,1];
tic
[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');
time = toc
disp(' ');
disp('The Result of Bernstein Computation of exemple-2-5-2 with real number approx is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
res = max(abs(bmin),abs(bmax));
disp(' ');