function [ bmin ] = matrix_min( bern_coeffs, params )
    expr = bern_coeffs(1);
    [c,t] = coeffs(expr,params);
    bmin = sum(abs(c)); % bmin is the bernstein min
    if isa(t(1),'double') % case where there is a lonely integer
        if (t(1)<0)
            bmin = bmin+2*t(1);  
        end
    end
    current = 0;
    for k=2:length(bern_coeffs) % we go through bernstein coefficients
        expr = bern_coeffs(k);
        [c,t] = coeffs(expr,params); % c is the set of coefficient
                                     % we must take the abs sum
        current = sum(abs(c));
        if isa(t(1),'double') % case where there is a lonely integer
            if (t(1)<0)
                current = current+2*t(1);  
            end
        end
        if current < bmin
            bmin = current;
        end
        current = 0;
    end
end