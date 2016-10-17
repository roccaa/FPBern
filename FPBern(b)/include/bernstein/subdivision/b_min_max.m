function [ C ] = b_min_max( params, B, P )
%B_MIN_MAX find min and max of all control points on the parameter vetrices

    % params : list of parameters
    % B : matrix of control points
    % P : Parameter space

    %addpath('../classic/');

    [rP,~] = size(P);
    [rB,cB] = size(B);
    
    for i =1:rP
        args{i} = P(i,:);
    end
    
    %All the possible bounds combinations
    allbounds = allcomb(args{:});
    [rab,~] = size(allbounds);
    
    vals = zeros(1,rab);
    
    for i=1:rB
        for j=1:cB
            for k=1:rab
               vals(k) = subs(B(i,j),params,allbounds(k,:));
            end
            minV = min(vals);
            maxV = max(vals);
            C{i,j} = [minV maxV];
        end
    end
end

