function [ pi,params ] = parser_poly( poly, vars, params_table)
    [c,t] = coeffs(poly,vars);
    pi = 0;
    params = [];
    pos = 0;
    for k=1:length(c)
        % We get the degree of each monom
        degree = feval(symengine, 'degree', t(k));
        mult = 1;
        for i=1:(degree-1)
            mult = mult*(1+params_table(i+pos));
            params = [params [params_table(i+pos)]];
        end
        if not(pi == 0)
            pi = (pi+c(k)*t(k)*mult)*(1+params_table(pos+(degree)));
            params = [params [params_table(pos+(degree))]];
            pos = pos+degree;
        else
            pi = c(k)*t(k)*mult;
            pos = pos+degree-1;
        end
    end
end
