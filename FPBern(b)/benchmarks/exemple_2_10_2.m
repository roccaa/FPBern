% Time

% Declare variables, parameters, and polynomial
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16;
pi = (((((((((( 30 * x3 + ( 30 * x5 + ( 3 * x6 + ( 3 * x7 + ( 3 * x8 + ( 3 * x9 + ( 3 * x10 + ( 3 * x11 + ( 3 * x12 + ( 3 * x13 + ( 3 * x14 + ( 2 * x15 + ( 3 * x16))))))))))))) * x1 + (( 270 * x3 + ( 30 * x4 + ( 300 * x5 + ( 30 * x6 + ( 30 * x7 + ( 30 * x8 + ( 30 * x9 + ( 30 * x10 + ( 30 * x11 + ( 30 * x12 + ( 30 * x13 + ( 30 * x14 + ( 20 * x15 + ( 30 * x16)))))))))))))) * x2)) * x1 + (( 1080 * x3 + ( 270 * x4 + ( 1350 * x5 + ( 135 * x6 + ( 135 * x7 + ( 135 * x8 + ( 135 * x9 + ( 135 * x10 + ( 135 * x11 + ( 135 * x12 + ( 135 * x13 + ( 135 * x14 + ( 90 * x15 + ( 135 * x16)))))))))))))) * x2^2)) * x1 + (( 2520 * x3 + ( 1080 * x4 + ( 3600 * x5 + ( 360 * x6 + ( 360 * x7 + ( 360 * x8 + ( 360 * x9 + ( 360 * x10 + ( 360 * x11 + ( 360 * x12 + ( 360 * x13 + ( 360 * x14 + ( 240 * x15 + ( 360 * x16)))))))))))))) * x2^3)) * x1 + (( 3780 * x3 + ( 2520 * x4 + ( 6300 * x5 + ( 630 * x6 + ( 630 * x7 + ( 630 * x8 + ( 630 * x9 + ( 630 * x10 + ( 630 * x11 + ( 630 * x12 + ( 630 * x13 + ( 630 * x14 + ( 420 * x15 + ( 630 * x16)))))))))))))) * x2^4)) * x1 + (( 3780 * x3 + ( 3780 * x4 + ( 7560 * x5 + ( 756 * x6 + ( 756 * x7 + ( 756 * x8 + ( 756 * x9 + ( 756 * x10 + ( 756 * x11 + ( 756 * x12 + ( 756 * x13 + ( 756 * x14 + ( 504 * x15 + ( 756 * x16)))))))))))))) * x2^5)) * x1 + (( 2520 * x3 + ( 3780 * x4 + ( 6300 * x5 + ( 630 * x6 + ( 630 * x7 + ( 630 * x8 + ( 630 * x9 + ( 630 * x10 + ( 630 * x11 + ( 630 * x12 + ( 630 * x13 + ( 630 * x14 + ( 420 * x15 + ( 630 * x16)))))))))))))) * x2^6)) * x1 + (( 1080 * x3 + ( 2520 * x4 + ( 3600 * x5 + ( 360 * x6 + ( 360 * x7 + ( 360 * x8 + ( 360 * x9 + ( 360 * x10 + ( 360 * x11 + ( 360 * x12 + ( 360 * x13 + ( 360 * x14 + ( 240 * x15 + ( 360 * x16)))))))))))))) * x2^7)) * x1 + (( 270 * x3 + ( 1080 * x4 + ( 1350 * x5 + ( 135 * x6 + ( 135 * x7 + ( 135 * x8 + ( 135 * x9 + ( 135 * x10 + ( 135 * x11 + ( 135 * x12 + ( 135 * x13 + ( 135 * x14 + ( 90 * x15 + ( 135 * x16)))))))))))))) * x2^8)) * x1 + (( 30 * x3 + ( 270 * x4 + ( 300 * x5 + ( 30 * x6 + ( 30 * x7 + ( 30 * x8 + ( 30 * x9 + ( 30 * x10 + ( 30 * x11 + ( 30 * x12 + ( 30 * x13 + ( 30 * x14 + ( 20 * x15 + ( 30 * x16)))))))))))))) * x2^9)) * x1 + (( 30 * x4 + ( 30 * x5 + ( 3 * x6 + ( 3 * x7 + ( 3 * x8 + ( 3 * x9 + ( 3 * x10 + ( 3 * x11 + ( 3 * x12 + ( 3 * x13 + ( 3 * x14 + ( 2 * x15 + ( 3 * x16))))))))))))) * x2^10);
params = [x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16];
vars = [x1 x2];
X = [-1,1,-1,1];
tic
[bmax,bmin] = bernstein_method(pi,params,vars,X,'polynomial');
time = toc
disp(' ');
disp('The Result of Bernstein Computation of exemple-2-10-2 with real number approx is :');
disp(' ');
disp('bmax=');disp(bmax);
disp('bmin=');disp(bmin);
res = max(abs(bmin),abs(bmax));
disp(' ');