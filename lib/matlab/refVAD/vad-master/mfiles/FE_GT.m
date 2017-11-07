function [gf,logE,X]=FE_GT(sam,fs,NbCh,FrameLen,tGtMat)

Efloor=exp(-50);
Gfloor=exp(-50);
PreEmp=0.97;

if ~exist('fs', 'var')
   fs=16000; % sampling frequency
end
if ~exist('FrameLen', 'var')
   FrameLen=0.025*fs;
end
FrameShift=0.01*fs;
Nfft=2^ceil(log(FrameLen)/log(2));

if ~exist('NbCh', 'var')
    NbCh=64; % number of GT filter banks
end
if ~exist('tGtMat', 'var')
    tGtMat=gammatone_matrix(Nfft,fs,NbCh)';
end
% truncate as in ReadWave
NbFr=floor( (length(sam)-FrameLen+FrameShift)/FrameShift);
sam=sam(1:NbFr*FrameShift+FrameLen-FrameShift);

% DC removal
sam=filter([1 -1],[1 -0.999],[0 reshape(sam,1,length(sam))]);

% framing
ind1=1:FrameShift:length(sam)-1-FrameLen+FrameShift;
ind2=(1:FrameLen)';
x=sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr)));

% logE
logE=log(max(sum(x.^2,1),Efloor));

% preemphasis & windowing
T=length(sam);
sam=[0 sam(2:T)-PreEmp*sam(1:T-1)];
win=hamming(FrameLen);
% Fourier Transform
X=fft(win(:,ones(1,NbFr)).*sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr))),Nfft);
X(1,:)=0;
X=abs(X(1:Nfft/2,:));
gf=max(tGtMat'*X,Gfloor).^(1/3);
