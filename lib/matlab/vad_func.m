function [ result, pp ] = vad_func( audio_dir, mode, threshold, output_type )
    
    disp('MRCG extraction ...')
    [data_len, winlen, winstep] = mrcg_extract( audio_dir );
    
    if mode == 3
        python_command = sprintf('python3 ./lib/python/VAD_test.py -m %d -l %d -b 100 --data_dir=./sample_data --model_dir=./saved_model --norm_dir=./norm_data', ... 
        mode, data_len);
    else
        python_command = sprintf('python3 ./lib/python/VAD_test.py -m %d -l %d -b 4096 --data_dir=./sample_data --model_dir=./saved_model --norm_dir=./norm_data', ... 
        mode, data_len);
    end
    mkdir './result'

    system(python_command);
    load('./result/pred.mat');
    
    pp = pred;
    result = zeros(length(pp), 1);
    result(pp>threshold) = 1;
    
    if output_type == 1
        result = frame2rawlabel(result, winlen, winstep);
        pp = frame2inpt(pp, winlen, winstep);
    end
     
end

