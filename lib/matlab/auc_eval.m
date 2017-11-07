function [ auc, pred, framed_label ] = auc_eval(label_dir)

load('./result/pred.mat');
load(label_dir);
audio_sr = 16000;
winlen             = ceil(audio_sr*25*0.001);	%window length (default : 25 ms)
winstep            = ceil(audio_sr*10*0.001);	%window step (default : 10 ms)
framed_label = Truelabel2Trueframe( y_label, winlen, winstep );
framed_label = framed_label(1:length(pred)); 
[~,~, ~, auc] = perfcurve(framed_label, pred, 1);
end

