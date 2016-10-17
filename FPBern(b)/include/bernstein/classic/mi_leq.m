function [ res ] = mi_leq( mi_a, mi_b )
%LEQ Check whether the multi-index mi_a is less equal than mi_b

    [~,na] = size(mi_a);
    [~,nb] = size(mi_b);
    
    if(na ~= nb)
        error('Multi-indexes have different sizes');
    end
    
    res = all(mi_a <= mi_b);
end

