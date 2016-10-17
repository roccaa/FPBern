function [ B ] = mat_contr_pts_dom( mi_lst, ap, X, options )
%MAT_CONTR_PTS compute Bernstein control points using a matrix approach
    % and mapping the domain

    [A,dims] = num2mat(mi_lst,ap);

    %addpath('../../multimatrix');

    [~,cA] = size(A);
    [~,cd] = size(dims);
    
    % Univariate polynomial
    if(cA == 1)
        B = Utilde(dims(1))*Vtilde(X(1,:),dims(1))*Wtilde(X(1,:),dims(1))*A;
        return;
    end
    
    
    B = A;
    
    for i=1:cd
        Ut = sym(Utilde(dims(1)));
        Vt = Vtilde(X(i,:),dims(1));
        Wt = Wtilde(X(i,:),dims(1));
        M = Ut*Vt*Wt;
        B = transp_2mat(M*B,dims+1);
        dims = [dims(2:end) dims(1)];
    end
    
    if strcmp(options,'m')
        return;
    else
        % Put all them in a row
        B = transpose(B);
        B = transpose(B(:));
    end
    
end

