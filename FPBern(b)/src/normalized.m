function [normalized_poly] = normalized(poly,vars,X)
%Normalize the polynome poly on [0,1]^m
k = 1; 
a = 0;
b = 0;
j = 0;
normalized_poly = poly;
    while (k < length(X)) 
        a = X(k);
        b = X(k+1);
        j = (k+1)/2;
        normalized_poly = subs(normalized_poly,vars(j),vars(j)*(b-a)+a);
        k = k+2;
    end
end

