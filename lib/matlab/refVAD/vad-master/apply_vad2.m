% ------------------------------------------------------------------------
% Copyright (C) 2016 University of Southern California, SAIL, U.S.
% Author: Maarten Van Segbroeck
% Mail: mvansegbroeck@gmail.com
% Date: 2016-20-12
% ------------------------------------------------------------------------
% vadout = apply_vad(audiodir,p1,p2);
%
%      Apply Voice Activity Detection to all files in a specified audio directory
%
%      --IN--
%      audiodir: directory of audio files (WAV format)
%      p1: speech/non-speech threshold [default:0.1]
%      p2: speech region smoothing [default:20]
%
%      --OUT--
%      vadout: VAD labels at frame level (10 ms)
%
%      When using please cite:
%
%      @article{van2013robust,
%       title={A Robust Frontend for VAD: Exploiting Contextual, Discriminative and Spectral Cues of Human Voice},
%       author={Van Segbroeck, Maarten and Tsiartas, Andreas and Narayanan, Shrikanth},
%       year={2013}
%     }

function vadout=apply_vad2(sam,fs_orig,p1,p2)
% 
% if ~exist('audiodir','var')
%   fprintf('Please specify an audio directory\n')
%   return
% end

% ---vad parameters (to tune)---
% speech noise threshold
if ~exist('p1','var')
  p1=0.003;  % default : 0.1
end
% speech region smoothing
if ~exist('p2','var')
  p2=60;  % default : 20
end

% set path
% addpath mfiles/


fs=8000;

% ---feature---
% GT
NbCh=64;
% Gabor
nb_mod_freq=2;
% LTSV
R=50; % context 50
M=10; % smoothing
ltsvThr=0.5;
ltsvSlope=0.2;
% vprob2 and ltsv2
K=30; order=4; % 30, 4
% -------------

% ---vad model---
load('models/model.mat')


% ---visualize---
visualize=false;



 % read in audio

 sam_8k=downsample(sam(:,1),fs_orig/fs);

 % [1] extract cochleagram
 gt=FE_GT(sam_8k,fs,NbCh);

 % [2] Gabor filtering applied on mel
 gbf=FE_GBF(sam_8k,fs,nb_mod_freq,false);
 gbf= [gbf gbf(:,ones(1,10)*size(gbf,2))];
 gbf = gbf(:,1:size(gt,2));

 % [3] LTSV
 ltsv=FE_LTSV(sam_8k,fs,R,M,gt,ltsvThr,ltsvSlope);
 ltsv2 = convert_to_context_stream(ltsv, K, order);
 ltsv2= [ltsv2 ltsv2(:,ones(1,10)*size(ltsv,2))];
 ltsv2 = ltsv2(:,1:size(gt,2));

 % [4] vprob prob
 vprob=voicingfeature(sam_8k,fs);
 vprob2 = convert_to_context_stream(vprob, K, order);
 vprob2 = [vprob2 vprob2(:,ones(1,10)*size(vprob,2))];
 vprob2 = vprob2(:,1:size(gt,2));

 % feature for VAD
 test_x = [gt;gbf;ltsv2;vprob2];
 test_x_norm = mvn(test_x);

 % VAD decoding
 [~,~,output] = nntest(dnn, test_x_norm');
 outprob=double(output(:,1));
 vadout=medfilt1(outprob.^2,p2)>p1;

 if visualize
  imagesc(mvn(gt));axis xy;hold on;
  plot(10*vadout,'m','LineWidth',3); zoom xon; hold off
 end

