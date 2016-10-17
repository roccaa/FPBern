dims = [7 7 7 7 7 7 7];
M = rand(size(alloc_mat(dims)));

tic;
navie_transp_2mat(M,dims);
fprintf('naive: '); 
toc;

tic;
transp_2mat(M,dims);
fprintf('dre: '); 
toc;