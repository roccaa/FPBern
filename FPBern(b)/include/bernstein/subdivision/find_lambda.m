function [ lambda ] = find_lambda(B, degs, P, params, r )
%FIND_LAMBDA

    [rB,cB] = size(B);

    C = b_min_max(params,B,P);
    
    X = [];
    Y = [];
    for i=1:rB
        for j=1:cB
            b = [i j] - 1;
            a = t2n(b,degs+1);
            ix = a(r);
            minMax = C{i,j};            
            X = [X ix];
            Y = [Y minMax(1,1)];
            X = [X ix];
            Y = [Y minMax(1,2)];            
        end
    end
    
    if(all(Y) > 0)
        error('All the control points are positive. Subdivision is useless.');
    end
    
    % Scale X to the unit box
    X = X/degs(r);
    
    % Polytope coordinates
    Poly = unique([X' Y'],'rows');
    % Convex hull of Poly
    K = convhull(Poly(:,1),Poly(:,2));
    % and its coordinates
    H = Poly(K,:);
    
    % Find the intersection with y-axis
    [rH,~] = size(H);
    Z = [];
    for i=1:rH-1
        if( (H(i,2) > 0 && H(i+1,2) <= 0 ) || (H(i+1,2) > 0 && H(i,2) <= 0 ))
            Z = [Z ; H(i,:) ; H(i+1,:)];
        end       
    end
    
    % Find intersections with Y-axis
    mu1 = polyxpoly([Z(1,1) Z(2,1)],[Z(1,2) Z(2,2)],[0 1],[0 0]);
    mu2 = polyxpoly([Z(3,1) Z(4,1)],[Z(3,2) Z(4,2)],[0 1],[0 0]);
    
    if (1-mu2 > mu1)
        lambda = mu2;
    else
        lambda = mu1;
    end    
end

