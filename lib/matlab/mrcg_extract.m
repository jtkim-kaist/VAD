function [data_len, winlen, winstep] = mrcg_extract( audio_dir, varargin )
    
    %UNTITLED2 이 함수의 요약 설명 위치
    %   자세한 설명 위치
    [noisy_speech, audio_sr] = audioread(audio_dir);
    
    if isempty(varargin)
        y_label = zeros(length(noisy_speech), 1);
    else
        label_dir = varargin{1};
        load(label_dir);
    end
    
    
    mkdir './sample_data'
    mkdir './sample_data/Labels'

    save_dir = './sample_data';
    name_mrcg = [save_dir, '/mrcg'];
    name_label = [save_dir, '/Labels/label'];
    
    
%     noisy_speech = noisy_speech(1:1000000,:);
%     y_label = y_label(1:1000000,:);
    
    mrcg_mat = MRCG_features(noisy_speech, audio_sr)';
    winlen             = ceil(audio_sr*25*0.001);	%window length (default : 25 ms)
    winstep            = ceil(audio_sr*10*0.001);	%window step (default : 10 ms)
    train_mean = mean(mrcg_mat,1);
    train_std = std(mrcg_mat,1,1);
    
    %% Read label
    
    framed_label = Truelabel2Trueframe( y_label, winlen, winstep );
    num = 0;
    
    if (length(mrcg_mat) > length(framed_label))
        binary_saver( name_mrcg, mrcg_mat(1:length(framed_label), :), num );
        binary_saver( name_label, framed_label, num );
        data_len = length(framed_label);
    else
        binary_saver( name_mrcg, mrcg_mat, num );
        binary_saver( name_label, framed_label(1:length(mrcg_mat), 1), num );
        data_len = length(mrcg_mat);
    end
    
    %% Save global normalization factor

    save([save_dir, '/normalize_factor'], 'train_mean', 'train_std');
    
    disp('MRCG extraction is successifully done.')
    
    
end

