function [bmax,bmin] = bernstein_method(expr,params,vars,X,mode)
    if (strcmp(mode,'polynomial'))
        disp(' ');
        disp('****************************');
        disp('***** POLYNOMIAL MODE ******');
        disp('****************************');
        disp(' ');
        disp('Load...');
        [bmax,bmin] = bernstein_method_polynomial(expr,params,vars,X);
    elseif (strcmp(mode,'rational'))
        disp(' ');
        disp('****************************');
        disp('****** RATIONAL MODE *******');
        disp('****************************');
        disp(' ');
        disp('Load...');
        [p,q] = numden(expr);
        [bmax,bmin] = bernstein_method_rational(p,q,params,vars,X);
    elseif (strcmp(mode,'example'))
        disp(' ');
        disp('****************************');
        disp('******* EXAMPLE MODE *******');
        disp('****************************');
        disp(' ');
        disp('We can for example use the bernstein method on f : (x1,x2) -> x1+x1*x2*p in [-1,1]*[-2,2], p is a parameter');
        disp('In that case we enter : bernstein_method(x1+x1*x2*p,[p],[x1 x2],[-1,1,-2,2],"polynomial")');
        disp('It will return the result of the computation and the result without the multiplication per epsilon.');
        disp(' ');
        disp('A Lot of example are available :');
        disp('sqroot');
        disp('sineTaylor');
        disp('sineOrder3');
        disp('rigibody1|rigibody2');
        disp('kepler0|kepler1|kepler2');
        disp('himmilbeau');
        disp('rational_g|rational_f');
        disp(' ');
        disp('To run one of them just use example where example is the name of the example above');
        disp('To run all use bern');
        disp(' ');
    else
        disp(' ');
        disp('Enter a correct mode : polynomial | rational | example');
        disp(' ');
        disp('Usage : bernstein_method(<expr>,<params>,<vars>,<X>,<mode>)');
        disp(' ');
    end
end