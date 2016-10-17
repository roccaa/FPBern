function [ isSharp ] = check_sharpness( b_contr_pts, P, params )

        isSharp = 0;
        C = b_min_max(params,b_contr_pts,P);
        [rC,cC] = size(C);
        maxB = -Inf;
        for i=1:rC
            for j=1:cC
                cp = C{i,j};
                cp(1,2);
                maxB = max(maxB,cp(1,2));
            end
        end
        
        cp = C{1,1};
        if ( maxB == cp(1,2) &&  cp(1,2) > 0)
            isSharp = 1;
        end
        cp = C{rC,cC};
        if ( maxB == cp(1,2) &&  cp(1,2) > 0)
            isSharp = 1;
        end
end

