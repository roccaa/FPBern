% Time

% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11;
pi = (( 12 * x3 + ( 12 * x5 + ( 6 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11)))))))) * x1 + (( 12 * x3 + ( 12 * x4 + ( 24 * x5 + ( 12 * x6 + ( 4 * x7 + ( 6 * x8 + ( 8 * x9 + ( 10 * x10 + ( 12 * x11))))))))) * x2)) * x1 + (( 12 * x4 + ( 12 * x5 + ( 6 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11)))))))) * x2^2);
params = [x3 x4 x5 x6 x7 x8 x9 x10 x11];
vars = [x1 x2];
X = [-1,1,-1,1];
tic
[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');
time = toc
disp(' ');
disp('The Result of Bernstein Computation of exemple-2-2-5 with real number approx is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
res = max(abs(bmin),abs(bmax));
disp(' ');