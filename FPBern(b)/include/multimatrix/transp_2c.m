function [ bt ] = transp_2c( b, d )
%TRANSPOSE a 2 dimenstional coordinate

    [~,cd] = size(d);

    bt(1,1) = mod(b(1,2),d(1,2));
    bt(1,2) = (b(1,2)-(mod(b(1,2),d(1,2))))/d(1,2) + b(1,1)*prod(d(3:cd));
    
end

