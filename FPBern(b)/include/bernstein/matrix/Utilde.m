function [ Ut ] = Utilde( n )
%UTILDE compute the inverse of U matrix

    Ut(1:n+1,1) = 1;
    Ut(1,2:n+1) = 0;
    Ut(n+1,2:n+1) = 1;
    
    for i=1:n-1
        for j=1:i
            Ut(i+1,j+1) = nchoosek(i,i-j)/nchoosek(n,j);
            Ut(i+1,i+2:n+1) = 0;
        end
    end


end

