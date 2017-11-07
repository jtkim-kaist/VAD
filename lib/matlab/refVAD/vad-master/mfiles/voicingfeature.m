%------------------------------------------------------------------------
% Copyright (C) 2012 University of Southern California, SAIL, U.S.
% Author: Maarten Van Segbroeck
% Mail: maarten@sipi.usc.edu
% ------------------------------------------------------------------------

%% DENOISING OF SPEECH
% i. removing of the voicing information contained in the signal
% ii. estimating the noise by minimum statistics of signal with removed voicing
% iii. voicing controlled spectral subtraction of noise from the original signal
%
function [vuv,mbvuv,correlogram]=voicingfeature(sam,fs,vThr,vSlope,verbose,dn,p,po)
% --IN--
% sam: signal
% sam_aper: voicing removed signal (optional and for speed up, use [] if unavailable)
% fs: sampling frequency
% p: noise tracking parameters (described below)
% po: spectral subtraction (described below)
% --OUT--
% dnsam: denoised signal
% sam_aper: voicing removed signal

% Show progress
if ~exist('verbose', 'var'),
    verbose=false;
end
if ~exist('dn', 'var'),
    dn=false;
end
if ~exist('vThr', 'var'),
    vThr=0.5;
end
if ~exist('vSlope', 'var'),
    vSlope=0.1;
end

% Default settings for noise tracking
if ~exist('p', 'var'),
    p = [4 0.2 1.0 1.0 1.0 1.0 0];
end
vsm = p(1); % size of running window for minimum statistics noise estimation (20)
nfl = p(2); % noise floor level of aperiodic signal (0.4): lower for low SNR
vln = p(3); % upper noise level of voicing (1.25)
vls = p(4); % lower speech level of voicing (1.5)
rfn = p(5); % noise reduction factor during noise frames (2.0)
rfs = p(6); % noise reduction factor during speech frames (0.5)
nsa = p(7); % add noise to compensate for introduced non-linear effect (0)

% Default settings for spectral subtraction
if ~exist('po', 'var'),
    po=[0.04 0.1 0.02 0.01 0.08 800 4 1.5 0.02 4]';
end
tg=po(1); % smoothing time constant for signal power estimate (0.04): high=reverberant, low=musical
ta=po(2); % smoothing time constant for signal power estimate (0.1)
tw=po(3); % frame/fft window length in ms (will be rounded up to a power of 2 samples) (0.025)
ts=po(4); % frame shift in ms (0.01)
to=po(5); % time constant for oversubtraction factor (0.08)
fo=po(6); % oversubtraction corner frequency (800)
km=po(7); % number of minimisation buffers to use (4)
kn=po(8); % noise estimate compensation (1.5)
kf=po(9); % subtraction floor (0.02)
ko=po(10); % oversubtraction scale factor (4)

ns=length(sam);

nw=pow2(nextpow2(fs*tw/km))*km;
nw=tw*fs;
ni=fs*ts;
ti=ni/fs;
nf=floor((ns-tw*fs+ni)/ni);

zg=exp(-ti/tg);
za=exp(-ti/ta);
zo=exp(-ti/to);

px=zeros(1+nw/2,1);
pxn=px;
os=px;
mb=ones(1+nw/2,km)*nw/2;
osf=ko*(1+(0:nw/2).'*fs/(nw*fo)).^(-1);

% remove voicing
[p,vuv,mbvuv,correlogram]=remove_voicing(sam,fs,tw,ts,verbose,vThr,vSlope);
return

%% REMOIVE VOICING
% i. estimate the pitch
% ii. pitch synchronous decomposition of the signal in its pitch harmonic component and aperiodic component
%
function [p,vuv,mbvuv,correlogram]=remove_voicing(sam,fs,tw,ts,verbose,vThr,vSlope)
% --IN--
% x: signal
% fs: sampling frequency
% tw: frame window length
% ts: frame window shift
% verbose: show progress
% --OUT--
% sam_aper: aperiodic part of the signal (voicing removed)

if ~exist('fs', 'var')
    fs = 16000;
end
if ~exist('tw', 'var'),
    tw = 0.025;
end
if ~exist('ts', 'var'),
    ts = 0.01;
end
if ~exist('verbose', 'var'),
    verbose = true;
end

% DC removal
sam=filter([1 -1],[1 -0.999],reshape(sam,1,length(sam)));

% preemphasis
PreEmp=0.97;
T=length(sam);
sam=[0 sam(2:T)-PreEmp*sam(1:T-1)];
% windowing
FrameLen=round(tw*fs);
FrameShift=round(ts*fs);
NbFr=floor( (T-FrameLen+FrameShift)/FrameShift);
ind1=(1:FrameShift:NbFr*FrameShift)-1;
ind2=(1:FrameLen)';
win=hamming(FrameLen);
xWin=single(win(:,ones(1,NbFr)).*sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr))));
% show progress
if verbose, fprintf('[1/3] pitch estimation .. '); end;

% pitch estimation
[p,vuv,mbvuv,correlogram]=pitch_estimation(xWin,fs,vThr,vSlope);
return

%% COMPUTE PITCH
% i. pitch estimation using a subharmonic decomposition method
%
function [pMicro,vuv,mbvuv,rfft,p,praw,score]=pitch_estimation(x,fs,vThr,vSlope,dpWidth)
% --IN--
% x: windowed and framed signal. One frame per column.
% fs: sampling frequency
% dpWidth: number of subharmonic bins (48 per octave) the pitch can change per frame
% --OUT--
% pMicro: pitch in Hz after micro adjustment and dp smoothing
% p: pitch after dp smoothing
% praw: subharmonic argmax; contains doubling/halving errors
% score: local score along optimal path in dp
% vuv: voiced/unvoiced measure. Beware: omitting this return value saves time

Nfft=1024;
if nargin<2,
  fs = 8000;
end
if ~exist('dpWidth', 'var'),
   dpWidth=4;
end

[FrameLen,T]=size(x);
interval = fs/Nfft;

XWin = fft(x,Nfft);
% keep only lower part of xfft
Nbw=floor(1250/fs*Nfft);
xfft = abs(XWin(1:Nbw,:));
flin=0:interval:interval*(Nbw-1); % frequency axis of xfft

% Find Relative Maximum
xAug=[zeros(1,T);xfft;zeros(1,T)];
maxpos = find( (xAug(2:end-1,:) > xAug(1:end-2,:)) & (xAug(2:end-1,:) > xAug(3:end,:)) );

% Set max positions and points around to zero
maxBar = xfft;
maxBar(maxpos)=0;
i=maxpos;i(rem(maxpos-1,Nbw)==0)=[];
maxBar(i-1)=0;
i=maxpos;i(rem(maxpos-1,Nbw)==Nbw-1)=[];
maxBar(i+1)=0;

% Subtract maxBar from x128 to get relative maxima and surounding points
relmax = xfft - maxBar;

% Low Pass Filter yeilding smoothed spectrum
smoothed = filter([1/4 1/2 1/4],1,[relmax zeros(Nbw,1)],[],2);
smoothed = smoothed(:,2:end);

% Spline Interpolation 48 points per octave
fsub = logspace(log10(2 * interval),log10(interval*Nfft),9*48+1);
fsub = fsub(fsub<=flin(Nbw) & fsub>=30);

% linear interpolation
iL=floor(fsub/interval)+1; % the lower index
wH=(fsub-flin(iL))'/interval;
wH=wH(:,ones(1,T));
subhar=(1-wH).*smoothed(iL,:)+wH.*smoothed(iL+1,:);

% Harmonic Summation
offset=round(48*log(1:15)/log(2)); % distance between harmonics
h=(0.84 .^ (0:14))';
[dum,Nmax]=min(abs(fsub-350)); % max pitch is 350 Hz
H=zeros(Nmax,T);
for s = 1:Nmax
    i = s + offset;
    i(i>size(subhar,1))=[];
    H(s,:)=sum(h(1:length(i),ones(1,T)).*subhar(i,:),1);
end
[dum,s_pitch]=max(H,[],1);
praw = fsub(s_pitch);
A=sparse(0.1*diag(ones(Nmax,1),0));
for k=-dpWidth:dpWidth,
    A=A+sparse(0.1*diag(ones(Nmax-abs(k),1),k));
end;
logA=spfun('log',A);
wght=2000;
[dum,s_pitch]=vit_gen(wght*H/max(max(H)),logA,1:Nmax,-1e-6*ones(Nmax,1),-wght/2);
score=H(sub2ind(size(H),s_pitch,1:T));
p = fsub(s_pitch);

% refine pitch estimate
pMicro=adjust_pitch(p/fs,XWin(2:Nfft/2,:),FrameLen)*fs;

if nargout>=2,
    Nover=2;
    rfft = real(ifft(xfft.^2,Nover*Nfft));
    Tsamples = 1./p .* fs;
    i = sub2ind(size(rfft),round(Nover*Tsamples)+1,1:T);
    vuv = rfft(i)./rfft(1,:);
    vuv=sigmoid(vuv-vThr,vSlope);
    mbvuv=[];
    Nover=1;
    NbMB=7;
    NbMB=1;
    for k=1:NbMB,
      i = sub2ind(size(rfft),round(k*Nover*Tsamples)+1,1:T);
      mbvuv(k,:) = rfft(i)./rfft(1,:);
    end

    %tGtMat=gammatone_matrix(Nfft,fs,Nfft/2)';
    %subplot 211; imagesc(rfft(1:Nfft/2,500:800),[-0.005 0.005]);axis xy
    %subplot 212; plot(tGtMat*rfft(1,500:800)); hold on; plot(smooth(rfft(i(500:800))./rfft(1,500:800),5),'r');hold off;axis tight
end

%% MODIFIED VITERBI
function [cost,path,fi,nActive]=vit_gen(logB,logA,PdfNr,fiInit,beam)
% --IN--
% logB: local scores for all Pdfs; size(logB,2)=T
% logA: transition matrix with log-probs (sparse)
% PdfNr: for every state, the index into logB to obtain state scores
% fiInit: initial state cost, sparse
% beam: (optional) normally negative
% --OUT--
% cost: total cost (scalar)
% path: stateNr (1-by-T)
% fi: state costs
% nActive: nnr active states

fi=sparse(fiInit);
N=length(logA);
[from,to,logAval]=find(logA);

T=size(logB,2);
Best=zeros(T,N);
nActive=N*ones(1,T);
for t=1:T,
   [fi,Best(t,:)]=spmax2(from,to,fi(from),logB(PdfNr(to),t)+logAval,N);
   if nargin>=5,
      sel=fi>max(fi)+beam;
      nActive(t)=sum(sel);
      fi(~sel)=-inf;
   end
end

% backtrace
[cost,path(1,T)]=max(fi);
cost=full(cost);
for t=T-1:-1:1,
   path(t)=Best(t+1,path(t+1));
end

%% ADJUST NORMALIZED PITCH ESTIMATE
function y=adjust_pitch(pNorm,XWin,Nwin)
% --IN--
% pNorm: normalized pitch
% XWin: fft spectrum, bins 2:Nfft/2
% Nwin: window lenght
% --OUT--
% y: adjusted normalized pitch

[N,T]=size(XWin);
Nfft=2*(N+1);
R=ceil(2*Nfft/Nwin); % Fourier resolution in bins
rng=(-R:R)';
y=pNorm;
for t=1:T,
   bin=pNorm(t)*Nfft:pNorm(t)*Nfft:N-R;
   n=round(bin/pNorm(t)); % harmonic number
   nMax=max(n);
   offs=rng*(n/nMax);
   sel=round(offs+bin(ones(2*R+1,1),:));
   a=zeros(size(sel));
   a(:)=abs(XWin(sel,t));
   [dum,i]=max(sum(a,2));
   y(t)=y(t)+rng(i)/nMax;
end

%% VOICEBOX FILES
function y=rfft(x,n)
%RFFT     FFT of real data Y=(X,N)
% Data is truncated/padded to length N if specified.
%   N even:     (N+2)/2 points are returned with
%                       the first and last being real
%   N odd:      (N+1)/2 points are returned with the
%                       first being real



%      Copyright (C) Mike Brookes 1998
%
%      Last modified Fri Apr  3 14:57:19 1998
%
%   VOICEBOX home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
  y=fft(x);
else
  y=fft(x,n);
end
if size(y,1)==1
  m=length(y);
  y(floor((m+4)/2):m)=[];
else
  m=size(y,1);
  y(floor((m+4)/2):m,:)=[];
end

function x=irfft(y,n)
%IRFFT    Inverse fft of a conjugate symmetric spectrum X=(Y,N)
% Y contains FIX(1+N/2) complex samples from the spectrum: if argument N
% is specified then Y will be truncated or padded accordingly
% IMPORTANT: If N is odd, it MUST be specified explicitly.
%
% See also RFFT



%      Copyright (C) Mike Brookes 1998
%
%      Last modified Fri Apr  3 14:57:14 1998
%
%   VOICEBOX home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fl=size(y,1)==1;
if fl y=y(:); end
[m,k]=size(y);
if nargin<2 n=2*m-2;
else
  mm=1+fix(n/2);
  if mm>m y=[y; zeros(mm-m,k)];
  elseif mm<m y(mm+1:m,:)=[];
  end
  m=mm;
end
if rem(n,2)             % odd case
  x=real(ifft([y;conj(y(m:-1:2,:))]));
else                    % even case
  y(m,:)=real(y(m,:));  % force nyquist element real
  w=ones(1,k);
%  t=[cumprod([-0.5i; exp(2i*pi/n)*ones(m-2,1)]); 0.5i];
  t=-0.5i* exp((2i*pi/n)*(0:m-1)).';
  z=(t(:,w)+0.5).*(conj(flipud(y))-y)+y;
  z(m,:)=[];
  zz=ifft(z);
  x=zeros(n,k);
  x(1:2:n,:)=real(zz);
  x(2:2:n,:)=imag(zz);
end

if fl x=x.'; end
