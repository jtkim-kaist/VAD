clc
clear
close all

addpath(genpath('./lib'));

%% parameter setting

audio_dir = './data/clean/clean_speech.wav';
mode = 0;           % 0 : ACAM3, 1 : bDNN, 2 : DNN, 3 : LSTM
threshold = 0.5;    % threshold for hard decision
output_type = 1;    % 0 : frame based prediction, 1: sample based prediction

%% prediction
% result : binary decision
% pp : posterior probability

[result, pp] = vad_func(audio_dir, mode, threshold, output_type);

%% plot (sample based)

label_dir = './data/clean/clean_label.mat'; % groud truh directory
load(label_dir);
s = audioread(audio_dir);

figure

p1 = plot(s);
hold on
p2 = plot(y_label*0.5, 'g--') ;
p3 = plot(result*0.3, 'r');
ylim([-0.3 0.7]);
legend([p2 p3],'ground truth', 'prediction')