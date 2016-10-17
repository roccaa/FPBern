function [bmax,bmin] = bernstein_method_rational(p,q,params,vars,X)
    % epsilon machine
    epsilon = 2^(-53);
    
    % Normalization of p and q
    normalized_p = normalized(p,vars,X);
    normalized_q = normalized(q,vars,X);

    % Extract exponents and coefficients in power basis
    [exps, coeffs] = sym2num(normalized_p,vars,params);
    % Compute Bernstein coefficients (with matrix method)
    bern_coeffs_p = mat_contr_pts(exps,coeffs,'r');
    disp('Bernstein Coeff Computation of p over.');
    disp('Load...');
    
    % Extract exponents and coefficients in power basis
    [exps, coeffs] = sym2num(normalized_q,vars,params);
    % Compute Bernstein coefficients (with matrix method)
    bern_coeffs_q = mat_contr_pts(exps,coeffs,'r');
    disp('Bernstein Coeff Computation of q over.');
    disp('Load...');
    
    if ~(isequaln(params,[]))
        [n m] = size(transpose(bern_coeffs_p));
        [n2 m2] = size(power(bern_coeffs_q,-1));
        if (m == n2) 
            matrix = transpose(bern_coeffs_p)*power(bern_coeffs_q,-1);
        elseif (n == m2)
            matrix = power(bern_coeffs_q,-1)*transpose(bern_coeffs_p);
        elseif (n == n2)
            matrix = bern_coeffs_p*power(bern_coeffs_q,-1);
        else
            matrix = power(bern_coeffs_q,-1)*bern_coeffs_p;
        end
        % Maximum in absolute of the bern_coeffs (with e in [-1;1])
        maximum = matrix_max(matrix, params);
        % Minimum in absolute of the bern_coeffs (with e in [-1;1])
        minimum = matrix_min(matrix, params);
        % Result
        bmax = double(epsilon*maximum);
        bmin = double(epsilon*minimum);
    else 
        % No params
        matrix = transpose(bern_coeffs_p)*power(bern_coeffs_q,-1);
        % Result
        bmax = double(max(max(matrix)));
        bmin = double(min(min(matrix)));
    end    
end