function [ exps, c  ] = sym2num( poly, vars, params )
%SYM2NUM Transform a symbolic polynomial into numeric representation

    [~,nvars] = size(vars);    
    [c,t] = coeffs(poly,vars);
    nter = size(t,2);
    exps = [];
    
    for i=1:nter
        if( t(i) ~= 1 )
            mono = sympoly(t(i) + t(i)*prod(vars));
            e = mono.Exponents(1,:);
            for j = 1:size(mono.Variables,2)
                exps(i,find(vars == mono.Variables{j})) = e(j);
            end
        else
            exps = [exps ; zeros(1,nvars)];
        end
    end
end

