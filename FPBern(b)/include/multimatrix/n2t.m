function [ b ] = n2t( a, d )
%N2T conver a n-dim coordinate to a 2-dim coordinate

    [~,ca] = size(a);
    [~,cd] = size(d);
    
    if(ca ~= cd)
        error('Coordinates and degree mismatching');
    end
    
    if(~all(a<=d))
        error('Coordinates cannot be lager than degrees');
    end
    
    b(1,1) = a(1,1);
    b(1,2) = a(1,2);
    
    for i=3:ca
        %prod(d(2:i-1))
        b(1,2) = b(1,2) + (a(1,i)*prod(d(2:i-1)));
    end
end

