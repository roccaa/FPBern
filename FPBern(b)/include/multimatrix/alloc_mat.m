function [ M ] = alloc_mat( dims )
%ALLOC_MAC allocate space for a 2dim matrix

    r = dims(1);
    c = prod(dims(2:end));
    
    M = zeros(r,c);

end

