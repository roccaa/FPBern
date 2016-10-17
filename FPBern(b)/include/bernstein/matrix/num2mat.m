function [ A, d ] = num2mat( mi_lst, ap )
%NUM2NMAT create 2d matrix with the numerical polynomial coefficients 
    % and find the degrees
    
   % addpath('../../multimatrix/');
    
    [ri,ci] = size(mi_lst);
    
    % Fetch the polynomial degree
    d = max(mi_lst,[],1);
    
    A = sym(alloc_mat(d+1));
    
    % univariate cases
    if ci == 1
        for i=1:ri
            A(mi_lst(i,:)+1,1) = ap(i);
        end
    else
        for i=1:ri
            b = n2t(mi_lst(i,:),d+1);
            A(b(1,1)+1,b(1,2)+1) = ap(i);
        end
    end
    
end

