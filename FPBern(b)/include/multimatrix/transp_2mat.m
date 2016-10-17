function [ Mt ] = transp_2mat( M, dims )
%TRANSP_2MAT transpose a 2-dim matrix

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
    
    Mt = sym(alloc_mat([dims(2:end) dims(1)]));    
    %Mt = alloc_mat([dims(2:end) dims(1)]);    
    
    for i=1:n
        b = [r(i) c(i)] - 1;
        bt = transp_2c(b,dims) + 1;
        Mt(bt(1),bt(2)) = M(r(i),c(i));
    end
    

end