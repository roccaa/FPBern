function [ l,h ] = l_h_function( pi, vars, params)
    l = 0;
    h = 0;
    [c,t] = coeffs(pi,vars);
    for k=1:length(c)
        coeff_l = 0;
        coeff_h = 0;
        [c2,t2] = coeffs(c(k),params);
        for i=1:length(c2)
            degree = feval(symengine, 'degree', t2(i));
            if (degree > 1)
                coeff_h = coeff_h+c2(i)*t2(i);
            end
            if (degree == 1)
                coeff_l = coeff_l+c2(i)*t2(i);
            end
        end
        l = l+coeff_l*t(k);
        h = h+coeff_h*t(k);
    end
end
