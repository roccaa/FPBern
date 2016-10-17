function [ BA, BB ] = subdivide( r, degs, lambda, B )
%SUBDIVIDE refine Bernstein control points by subdivision
    % r : direction
    % degs : variable degrees
    % X : current box space
    % lambda: subdivision point
    % B : control points in X
    
    [~,rX] = size(degs);
    
    if (r <= 0 || r > rX)
        error('Subdivision direction r must be larger than 0 and smaller then box space');
    end
    
    if (lambda < 0 || lambda > 1)
        error('Subdivision point Lambda must be between 0 and 1');
    end
    
    
    [rB,cB] = size(B);
    
    % Refine the coefficients in the sub-box A (left of lambda)
    BA = B;
    BB = B;
    
    for mi=1:rB
        for mj=1:cB
            
            % Compute the current multi-index
            i = t2n([mi mj]-1,degs+1);
            
            for k=1:degs(r)
                for j=1:rX
                    if j ~= r
                        for ij = 0:degs(j)
                            
                            % Refine on the A-side
                            if ij<k                               
                                BA(mi,mj) = BA(mi,mj);
                            else
                                ni = i;
                                ni(j) = ij - 1;
                                bi = n2t(ni,degs+1)+1;
                                BA(mi,mj) = (1-lambda)*BA(bi(1,1),bi(1,2)) + lambda*BA(mi,mj);
                            end
                            
                            % Refine on the B-side
                            ni = i;
                            ni(j) = degs(j) - k;
                            bi = n2t(ni,degs+1)+1;
                            BB(bi(1,1),bi(1,2)) = BA(mi,mj); 
                            
                        end
                    end
                end
            end           
            
        end
    end
    
end

