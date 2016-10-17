function [ res ] = contr_pts_degs( mi_lst,ap,degs,options)
% CONTR_PTS Multi-index Bernstein control points
    % mi_lst: multi-indexes list
    % ap: coefficient list
    % options: 'm' to have the coeffs in matrix form
    
    or_mi_lst = mi_lst;
    
    for i=1:prod(degs+1)
        val = i-1;
        e = [];
        for j=size(degs,2):-1:1        
            e(j) = mod(val,degs(j)+1);
            val = idivide(int16(val),int16(degs(j)+1));
        end
        if ~(is_in(e,or_mi_lst))
            mi_lst = [mi_lst ; e];
            ap = [ap 0];
        end
    end

    [~,cidx]=size(mi_lst); 
    mi_lst = [mi_lst ; zeros(1,cidx)];
    real_list = cell(cidx,1);
    
    
    % Fetch the polynomial degree
    d = max(mi_lst);
    
    % Allocate space for all the multi-index combinations
    for i=1:cidx
        real_list{i} = 0:1:d(1,i);
    end
    
    %tic;
    % Compute all the combinations
    Id = allcomb(real_list{:});
    [ridx,~]=size(Id);
    [~,cidx]=size(ap); 
    
    syms real_ap;
    for j = 1:ridx
        real_ap(1,j) = 0;
    end
    
    % Associate each multi-index with its coefficient
    for i=1:cidx
        indx = find(ismember(Id,mi_lst(i,:),'rows'),1);
        real_ap(1,indx) = ap(1,i);
    end
    %fprintf('   Create multi-index: ');toc;
    
    %tic;
    % Compute Bernstein coefficients 
    if(strcmp(options,'m'))
        addpath('../../multimatrix/');
        res = alloc_mat(d+1);
        for i=1:ridx
            b = n2t(Id(i,:),d+1)+1;
            res(b(1,1),b(1,2)) = bip(Id(i,:),real_ap,Id,d);
        end        
    else
        for i=1:ridx
            res(i) = bip(Id(i,:),real_ap,Id,d);
        end
    end
    
    %fprintf('   Compute Bern: ');toc;
    
end

function [ res ] = is_in( row, mi_lst)
    for i=1:size(mi_lst,1)
        if all(row == mi_lst(i,:))
            res = 1;
            return;
        end
    end
    res = 0;
end

