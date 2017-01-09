addpath(genpath('benchmarks/'))
addpath(genpath('include/'))
addpath(genpath('src/'))

programs_names = {'simple_test';'rigidBody1';'rigidBody2';'kepler0real';'kepler1';'kepler2';'sineTaylor';'sineOrder3';'sqroot';'himmilbeau';'schwefel';'magnetism';'caprasse';'exemple_2_2_5';'exemple_2_2_10';'exemple_2_2_15';'exemple_2_2_20';'exemple_2_5_2';'exemple_2_10_2';'exemple_5_2_2'};%'exemple_10_2_2'};
files_names = {'simple_test';'rigibody1real';'rigibody2real';'kepler0real';'kepler1real';'kepler2real';'sineTaylorreal';'sineOrder3real';'sqrootreal';'himmilbeaureal';'schwefel';'magnetism';'caprasse';'exemple_2_2_5';'exemple_2_2_10';'exemple_2_2_15';'exemple_2_2_20';'exemple_2_5_2';'exemple_2_10_2';'exemple_5_2_2'};%'exemple_10_2_2'};
%programs_names = {'exemple_2_2_5';'exemple_2_2_10';'exemple_2_2_15';'exemple_2_2_20'};%'exemple_10_2_2'};
%files_names = {'exemple_2_2_5';'exemple_2_2_10';'exemple_2_2_15';'exemple_2_2_20'};%'exemple_10_2_2'};

result = cell(3 ,size(files_names,1)+1);
result{1,1} = 'Name  ';
result{2,1} = 'Roundofferror  ';
result{3,1} = 'Total time  ';

n_repeat = 1;
for i=1:size(files_names,1)
    ttime = 0;
    for j=1:n_repeat
        eval(files_names{i});
        ttime = ttime + time;
    end    
    result{1,i+1} = programs_names{i};
    result{2,i+1} = res;
    result{3,i+1} = ttime/n_repeat;
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
