function [ a ] = t2n( b,d )
%T2N convert a 2-dim coordinate into a n-dim coordinate

    [~,cb] = size(b);
    if( cb ~= 2)
        error('b must have 2 coordinates')
    end
    
    [~,cd] = size(d);
    a(1,1) = b(1,1);
    for i=cd:-1:2
        div = prod(d(2:i-1));
        a(1,i) = floor(b(1,2)/div);
        b(1,2) = mod(b(1,2),div);
    end

end

