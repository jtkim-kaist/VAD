clc
clear

addpath(genpath('./lib'));


%% parameter setting
SE_list = getAllFiles('.\Sound_Effect');
r = 3; 
Fs = 16000;
save_dir = ['.\Sound_Effect_', num2str(16000), '\Sound_Effect_'];
%% paramter setting end



for i=1:1:length(SE_list)
    i
    x = audioread(SE_list{i});
    x = mean(x, 2);  % 2 channel -> 1 channel
    x = x/max(abs(x));
    x = decimate(x, r);
    num_file = sprintf('%5.5d',i);
    
    file_name = [save_dir, num_file, '.wav'];
    audiowrite(file_name, x, Fs)
        
end