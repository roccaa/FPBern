function [ V ] = Wtilde( I,n )
%WTILDE

    [~,c] = size(I);
    if( c ~= 2 )
        error('I must be an interval [a b].');
    end
    if( I(1) > I(2) )
        error('Lower bound must be smaller than the upper bound.');
    end

    V = diag(ones(n+1,1),0);
    
    for i=1:n+1
        for j=i+1:n+1
            V(i,j) = nchoosek(j-1,i-1)*I(1)^(j-i);
        end
    end
    
    


end

