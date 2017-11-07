% ------------------------------------------------------------------------
% Copyright (C) 2012 University of Southern California, SAIL, U.S.
% Author: Maarten Van Segbroeck
% Mail: maarten@sipi.usc.edu
% ------------------------------------------------------------------------

%% LOG-MEL SCALED SPECTRAL FEATURES
function [lp,logE]=FE_MEL(sam,fs,tMelMat)
% --IN--
% sam: signal
% fs: sampling frequency
% tMelMat: MEL matrix (optional)
% --OUT--
% lp: log-mel-power of sam
% logE: log-energy of sam

if ~exist('fs', 'var')
    fs = 16000;
end

FrameLen=floor(0.025*fs);
FrameShift=0.01*fs;
Efloor=exp(-50);
MELfloor=exp(-50);
PreEmp=0.97;
Nfft=2^ceil(log(FrameLen)/log(2));
NbCh=floor(hz2meldm(fs/2))-1; % number of MEL filter banks 
NbCh=24;
if ~exist('tMelMat', 'var')
    tMelMat=mel_matssp(fs,NbCh,Nfft)';
end

% truncate
NbFr=floor( (length(sam)-FrameLen+FrameShift)/FrameShift);
sam=sam(1:NbFr*FrameShift+FrameLen-FrameShift);

% DC removal
sam=filter([1 -1],[1 -0.999],[0 reshape(sam,1,length(sam))]);

% framing
ind1=(1:FrameShift:NbFr*FrameShift)-1;
ind2=(1:FrameLen)';
x=sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr)));

% logE
logE=log(max(sum(x.^2,1),Efloor));

% preemphasis & windowing
T=length(sam);
sam=[0 sam(2:T)-PreEmp*sam(1:T-1)];
win=hamming(FrameLen);
xNoWin=sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr)));
xWin=win(:,ones(1,NbFr)).*xNoWin;

% logMel
X=fft(xWin,Nfft);
X(1,:)=0;
Xmag=abs(X(1:Nfft/2+1,:));
lp=log(max(tMelMat*Xmag,MELfloor));

%% MEL-SCALED FILTERBANK
function [M,fCen]=mel_matssp(fs,NbCh,Nfft)
% --IN--
% fs: sampling frequency
% NbCh: number of MEL filters
% Nfft: number of FFT points 
% --OUT--
% M: Nfft/2+1-by-NbCh matrix of triangular weights
% fcen: center frequencies of triangular weights

fCen=meldm2hz(1:NbCh);
fLow=meldm2hz(0:NbCh-1);
fHigh=meldm2hz(2:NbCh+1);

flsamp = max(hz2spc(fLow,fs,Nfft/2+1),0);
i=find(spc2hz(flsamp,fs,Nfft/2+1)<fLow);
flsamp(i)=flsamp(i)+1;

fhsamp = min(hz2spc(fHigh,fs,Nfft/2+1),Nfft/2);
i=find(spc2hz(fhsamp,fs,Nfft/2+1)>fHigh);
fhsamp(i)=fhsamp(i)-1;

TotLen=fhsamp-fhsamp+1;
M=sparse([],[],[],Nfft/2+1,NbCh,sum(TotLen));
for k=1:NbCh,
  f=spc2hz(flsamp(k):fhsamp(k),fs,Nfft/2+1);
  w=zeros(size(f));
  sel=f<=fCen(k);
  w(sel)=max(0,1+(f(sel)-fCen(k))/(fCen(k)-fLow(k)));
  sel=f>fCen(k);
  w(sel)=max(0,(f(sel)-fHigh(k))/(fCen(k)-fHigh(k)));
  M(flsamp(k)+1:fhsamp(k)+1,k)=w';
end
