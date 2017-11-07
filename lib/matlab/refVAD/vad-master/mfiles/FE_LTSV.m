% ------------------------------------------------------------------------
% Copyright (C) 2013 University of Southern California, SAIL, U.S.
% Author: Maarten Van Segbroeck
% Mail: maarten@sipi.usc.edu
% Date: 2013-28-10
% ------------------------------------------------------------------------
% ltsv = FE_LTSV(sam,fs,R,M,ltsvThr,ltsvSlope);
%      Generates a Long-Term Spectral Variability (LTSV) stream of speech compute on the 
%      Gammatone filtered time-frequency representation of speech (different from original
%      LTSV stream in [1]). The LTSV is also processed by a sigmoid function to assure
%      a probability value between 0 and 1 on the current signal  
%
%      --IN--
%      sam: sample file
%      fs: sampling frequency
%      R: context window parameter [default:50 
%      M: smoothing parameter [default:M]
%      ltsvThr: threshold of sigmoid function [default:0.5]
%      ltsvSlope: slope of sigmoid function [default:0.2]
%
%      --OUT--
%      ltsv: LTSV stream of the signal

function ltsv=FE_LTSV(sam,fs,R,M,gt,ltsvThr,ltsvSlope);

FrameLen=0.032*fs;
FrameShift=0.01*fs;
Nfft=2^ceil(log(FrameLen)/log(2));
Nfft=128;
S=FE_GT(sam,fs,Nfft,FrameLen);

if ~exist('R', 'var')
 R=50; 
end
if ~exist('M', 'var')
 M=10; 
end
if ~exist('ltsvThr', 'var')
 ltsvThr=0.5;
end
if ~exist('ltsvSlope', 'var')
 ltsvSlope=0.2;
end

append_ndx=size(S,2):-1:max(size(S,2)-20,1);
S=[S S(:,append_ndx)];

% frequency smoothing
S_M=[];
for k=1:size(S,1)
	S_M(k,:)=(smooth(S(k,:),M,'moving'));
end
% normalized spectrogram
S_R=[];
for t=1:size(S,2)
	ndx=[ones(1,R/2-t) max(1,t-R/2+1):min(size(S,2),t+R/2) ones(1,R/2+t-size(S,2))];
%	ndx=[t:min(size(S,2),t+R) ones(1,R+t-size(S,2))];
	S_R(:,t)=S_M(:,t)./sum(S_M(:,ndx),2);
end

% entropy measure of normalized spectrogram over R consecutive frames, ending at current frame
E_R=-100*S_R.*log(100*S_R);
L=var(E_R);

L=[];
N=[1];
skip=0;
for n=1:length(N)
  for i=1:N(n)
        fndx=round([Nfft*(i-1)/(2*N(n))+1:Nfft*i/(2*N(n))]);
	fndx=unique(min(max(fndx,1),Nfft));
	L(i+skip,:)=var(E_R(fndx,:));
  end
  skip=skip+N(n);
end

S=S(:,1:end-length(append_ndx));
L=L(:,1:end-length(append_ndx));
if exist('ltsvThr', 'var') && exist('ltsvSlope', 'var')
 L=sigmoid(2*L-ltsvThr,ltsvSlope);
end
ltsv=L;
