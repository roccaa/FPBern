function [ b ] = bip( i,ap,mi_lst,d )
% BIP i-th parametrized Bernstein coefficient
%   i: multi-index
%   ap: parameter/coefficient list
%   mi_lst: multi-index list
%   d: polynomial degree (multi-index)

    [nidx,~] = size(mi_lst);
    b = 0;
    
    for j=1:nidx
        if mi_leq(mi_lst(j,:),i)
           b = b + (mi_nchoosek(i,mi_lst(j,:))/mi_nchoosek(d,mi_lst(j,:)))*ap(j);
        end        
    end

end

