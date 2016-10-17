function [ bmax ] = matrix_max( bern_coeffs, params )
    bmax = 0; % bmax is the bernstein max
    current = 0; % current is the current value 
    for k=1:length(bern_coeffs) % we go through bernstein coefficients
        expr = bern_coeffs(k);
        [c,t] = coeffs(expr,params); % c is the set of coefficient
                                     % we must take the abs sum
        current = sum(abs(c));
        if isa(t(1),'double') % case where there is a lonely integer
            if (t(1)<0)
                current = current+2*t(1);  
            end
        end
        if current > bmax
            bmax = current;
        end
        current = 0;
    end
end

