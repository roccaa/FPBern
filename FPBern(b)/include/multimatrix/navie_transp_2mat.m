function [ Mt ] = navie_transp_2mat( M,dims )
%NAVIE_TRANSP_2MAT Summary of this function goes here
%   Detailed explanation goes here
    
    [rM,cM] = size(M);
    r = [];
    c = [];
    
    % Fetch non-zero elements
    % we can't use [r,c] = find(M); if we want symbolic computations
    for i=1:rM
        for j=1:cM
            if (M(i,j) ~= 0)
                r = [r ; i];
                c = [c ; j];
            end
        end
    end

    [n,~] = size(r);
    
    %Mt = sym(alloc_mat([dims(2:end) dims(1)]));    
    Mt = alloc_mat([dims(2:end) dims(1)]);    
    
    for i=1:n
        
        b = [r(i) c(i)] - 1;
        a = t2n(b,dims);
        a = [a(2:end) a(1)];
        bt = n2t(a,[dims(2:end) dims(1)])+1;
        Mt(bt(1),bt(2)) = M(r(i),c(i));
    end

end

