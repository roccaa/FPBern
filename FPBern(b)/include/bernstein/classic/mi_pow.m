function [ res ] = mi_pow( B, midx )
% mi_exp multi-index power
%   B: base ([b1 b2 b3 ...])
%   midx: multi-index

    res = prod(power(B,midx));

end

