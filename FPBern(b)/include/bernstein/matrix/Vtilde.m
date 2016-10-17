function [ V ] = Vtilde( I, n )
%VTILDE

    [~,c] = size(I);
    if( c ~= 2 )
        error('I must be an interval [a b].');
    end
    if( I(1) > I(2) )
        error('Lower bound must be smaller than the upper bound.');
    end

    V = zeros(n+1,n+1);
    V(1,1) = 1;
    
    for i=2:n+1
        V(i,i) = (I(2) - I(1))^(i-1);
    end
end

