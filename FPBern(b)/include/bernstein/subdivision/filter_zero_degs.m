function [ B,C ] = filter_zero_degs( degs )
%FILTER_ZERO_DEGS

    [~,nd] = size(degs);
    C = zeros(1,nd);
    j = 1;
    for i=1:nd
        if(degs(1,i) ~= 0)
            B(1,j) = degs(1,i);
            C(1,i) = j;
            j = j+1;
        end
    end
end

