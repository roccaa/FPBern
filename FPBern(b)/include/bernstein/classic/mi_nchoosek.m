function [ res ] = mi_nchoosek( mi,md )
%MINCHOOSEK multi-index i choose d

    [~,ni] = size(mi);
    [~,nd] = size(md);
    
    if(ni ~= nd)
        error('Multi-indexes have different sizes');
    end
    
    res = 1;
    for i=1:ni
        res = res*nchoosek(mi(i),md(i));
    end


end

