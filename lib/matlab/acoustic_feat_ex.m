function acoustic_feat_ex( data_dir, save_dir )

rng(0);
%% Directory setting

% system(['rm -rf ', save_dir]);
% 
% system(['mkdir ', save_dir]);
system(['mkdir ', save_dir, '/Normalize_Factor']);
system(['mkdir ', save_dir, '/Labels']);

%% Parameter setting

audio_sr = 16000;
split_num = 1;
name_mrcg = [save_dir, '/se_mrcg'];
name_label = [save_dir, '/se_label'];

audio_list = getAllFiles(data_dir, 'FileFilter', '\.wav$');
label_list = getAllFiles(data_dir, 'FileFilter', '\.mat$');

winlen             = ceil(audio_sr*25*0.001);	%window length (default : 25 ms)
winstep            = ceil(audio_sr*10*0.001);	%window step (default : 10 ms)

train_mean = 0;
train_std = 0;

for i = 1:1:length(audio_list)
    clc
    fprintf("MRCG extraction %d/%d ...\n", i, length(audio_list));
    %% Read audio
    
    noisy_speech = audioread(audio_list{i});  % noisy_speech load
    noisy_speech = noisy_speech(1:(length(noisy_speech)-mod(length(noisy_speech), split_num)));
    noisy_speech = reshape(noisy_speech, [], split_num);
    
    %% Caliculate MRCG
    mrcg = cell(split_num, 1);
    
    for j = 1:1:split_num
        mrcg{j, 1} = MRCG_features(noisy_speech(:, j), audio_sr)';
        %     imagesc(s(20000:20500,:)*1000)
    end
    
    mrcg_mat = cell2mat(mrcg);
    
    size(mrcg_mat)
    %% Save normalization factor
    
    temp_mean = mean(mrcg_mat,1);
    temp_std = std(mrcg_mat,1,1);
    save([save_dir, '/Normalize_Factor/normalize_factor_', sprintf('%3.3d', i)],'temp_mean', 'temp_std');
    train_mean = temp_mean + train_mean;
    train_std = temp_std + train_std;
    
    %% Read label
    label = cell2mat(struct2cell(load(label_list{i})));  % label load
    
    %% Save framed label & MRCG
    framed_label = Truelabel2Trueframe( label, winlen, winstep );
    length(framed_label)
    if (length(mrcg_mat) > length(framed_label))
        binary_saver( name_mrcg, mrcg_mat(1:length(framed_label), :), i );
        binary_saver( name_label, framed_label, i );
    else
        binary_saver( name_mrcg, mrcg_mat, i );
        binary_saver( name_label, framed_label(1:length(mrcg_mat), 1), i );
    end
end

disp('MRCG extraction done.')
%% Save global normalization factor

global_mean = train_mean / length(audio_list);
global_std = train_std / length(audio_list);
save([save_dir, '/global_normalize_factor'], 'global_mean', 'global_std');

%% Move label data

feat_list = getAllFiles(save_dir);

for i=1:1:length(feat_list)
    if ~isempty(strfind(feat_list{i}, 'label'))
        [pathstr, name, ext] = fileparts(feat_list{i});
        new_path = [pathstr, '/Labels/', name, ext];
        movefile(feat_list{i}, new_path);
    end
end

end



