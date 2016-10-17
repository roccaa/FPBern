function [ bmax,bmin ] = matrix_max_min( bern_coeffs, params )
    bmax = 0; % bmax is the bernstein max
    bmin = 0; % bmax is the bernstein min
    expr = bern_coeffs(1);
    [c,t] = coeffs(expr,params);
    bmin = -sum(abs(c)); % bmin is the bernstein min
    bmax = -bmin;
    expr
    isempty(c)
    if isa(t(1),'double') % case where there is a lonely integer
        if (t(1)<0)
            bmax = bmax+2*t(1);
        else
            bmin = bmin+2*t(1);
        end
    end
    current = 0;
    imax = 1;
    for k=2:length(bern_coeffs) % we go through bernstein coefficients
        expr = bern_coeffs(k);
        [c,t] = coeffs(expr,params); % c is the set of coefficient
                                     % we must take the abs sum
        current_max = sum(abs(c));
        current_min = -current_max;
        if isa(t(1),'double') % case where there is a lonely integer
            if (t(1)<0)
                current_max = current_max+2*t(1);  
            else 
                current_min = current_min+2*t(1);
            end
        end
        if current_min < bmin
            imin = k;
            bmin = current_min;
        end
        if current_max > bmax
            imax = k;
            bmax = current_max;
        end
        current_max = 0;
        current_min = 0;
    end
    disp('imax = ');
    disp(imax);
end