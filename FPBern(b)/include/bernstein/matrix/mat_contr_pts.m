function [ B ] = mat_contr_pts( mi_lst, ap, options )
%MAT_CONTR_PTS compute Bernstein control points using a matrix approach

    [A,dims] = num2mat(mi_lst,ap);

    [~,cA] = size(A);
    [~,cd] = size(dims);
    
    % Univariate polynomial
    if(cA == 1)
        B = Utilde(dims(1))*A;
        return;
    end
    
    
    B = A;
    
    for i=1:cd
        Ut = sym(Utilde(dims(1)));
        B = transp_2mat(Ut*B,dims+1);
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

