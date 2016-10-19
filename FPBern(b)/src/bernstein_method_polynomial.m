function [bmax,bmin] = bernstein_method_polynomial(pi,params,vars,X)
    % epsilon machine   
    epsilon = 2^(-53); % 2^(-24) pour le flottant % previsously 2^(-53)

    % Normalization of pi
    normalized_pi = normalized(pi,vars,X);
    normalized_pi
    % Extract exponents and coefficients in power basis
    [exps, coeffs] = sym2num(normalized_pi,vars,params);

    % Compute Bernstein coefficients (with matrix method)
    bern_coeffs = mat_contr_pts(exps,coeffs,'r');
    disp('Bernstein Coeff Computation over.');
    disp('Load...');
    
    if ~(isequaln(params,[]))
        % Maximum in absolute of the bern_coeffs (with e in [-1;1])
        [bmax,bmin] = matrix_max_min(bern_coeffs, params);
        bmax = double(bmax*epsilon);
        bmin = double(bmin*epsilon);
    else
        bmax = max(bern_coeffs);
        bmin = min(bern_coeffs);
    end
end