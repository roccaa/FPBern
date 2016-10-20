addpath(genpath('benchmarks/'))
addpath(genpath('include/'))
addpath(genpath('src/'))

programs_names = {'simple_test';'rigidBody1';'rigidBody2';'kepler0real';'kepler1';'kepler2';'sineTaylor';'sineOrder3';'sqroot';'himmilbeau';'schwefel';'magnetism';'caprasse'};
files_names = {'simple_test';'rigibody1real';'rigibody2real';'kepler0real';'kepler1real';'kepler2real';'sineTaylorreal';'sineOrder3real';'sqrootreal';'himmilbeaureal';'schwefel';'magnetism';'caprasse'};

result = cell(3 ,size(files_names,1)+1);
result{1,1} = 'Name  ';
result{2,1} = 'Roundofferror  ';
result{3,1} = 'Total time  ';

for i=1:size(files_names,1)
    eval(files_names{i});
    result{1,i+1} = programs_names{i};
    result{2,i+1} = res;
    result{3,i+1} = time;
end

% simple_test
% rigibody1real
% rigibody2real
% kepler0real
% kepler1real
% kepler2real
% sineTaylorreal
% sineOrder3real
% sqrootreal
% himmilbeaureal
