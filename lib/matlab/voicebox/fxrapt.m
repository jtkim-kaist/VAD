function [fx,tt]=fxrapt(s,fs,mode,q)
%FXRAPT RAPT pitch tracker [FX,VUV]=(S,FS,M,Q)
%
% Input:   s(ns)      Speech signal
%          fs         Sample frequency (Hz)
%          mode       'g' will plot a graph [default if no output arguments]
%                     'u' will include unvoiced fames (with fx=NaN)
%          q          stucture with parameter values (e.g. q.f0min=40); see below for a list
%
% Outputs: fx(nframe)     Larynx frequency for each fram,e (or NaN for silent/unvoiced)
%          tt(nframe,3)  Start and end samples of each frame. tt(*,3)=1 at the start of each talk spurt
%
% Plots a graph if no outputs are specified showing lag candidates and selected path
%

% Bugs/Suggestions:
%   (1) Include backward DP pass and output the true cost for each candidate.
%   (2) Add an extra state to distinguish between voiceless and silent
%   (3) N-best DP to allow longer term penalties (e.g. for frequent pitch doubling/halving)

% The algorithm is taken from [1] with the following differences:
%
%      (a)  the factor AFACT which in the Talkin algorithm corresponds roughly
%           to the absolute level of harmonic noise in the correlation window. This value
%           is here calculated as the maximum of three figures:
%                   (i) an absolute floor set by p.absnoise
%                  (ii) a multiple of the peak signal set by p.signoise
%                 (iii) a multiple of the noise floor set by p.relnoise
%      (b) The LPC used in calculating the Itakura distance uses a Hamming window rather than
%          a Hanning window.
%
% A C implementation of this algorithm by Derek Lin and David Talkin is included as  "get_f0.c"
% in the esps.zip package available from http://www.speech.kth.se/esps/esps.zip under the BSD
% license.
%
% Refs:
%      [1]   D. Talkin, "A Robust Algorithm for Pitch Tracking (RAPT)"
%            in "Speech Coding & Synthesis", W B Kleijn, K K Paliwal eds,
%            Elsevier ISBN 0444821694, 1995

%      Copyright (C) Mike Brookes 2006-2013
%      Version: $Id: fxrapt.m 9312 2017-01-19 13:19:13Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
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
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=s(:); % force s to be a column
if nargin<4
    q=[];
    if nargin<3
        mode=' ';
    end
end
doback=0;   % don't do backwards DP for now

% set default parameters

p0.f0min=50;           % Min F0 (Hz)
p0.f0max=500;          % Max F0 (Hz)
p0.tframe=0.01;        % frame size (s)
p0.tlpw=0.005;         % low pass filter window size (s)
p0.tcorw=0.0075;       % correlation window size (s)
p0.candtr=0.3;         % minimum peak in NCCF
p0.lagwt=0.3;          % linear lag taper factor
p0.freqwt=0.02;        % cost factor for F0 change
p0.vtranc=0.005;       % fixed voice-state transition cost
p0.vtrac=0.5;          % delta amplitude modulated transition cost
p0.vtrsc=0.5;          % delta spectrum modulated transition cost
p0.vobias=0.0;         % bias to encourage voiced hypotheses
p0.doublec=0.35;       % cost of exact doubling or halving
p0.absnoise=0;         % absolute rms noise level
p0.relnoise=2;         % rms noise level relative to noise floor
p0.signoise=0.001;     % ratio of peak signal rms to noise floor (0.001 = 60dB)
p0.ncands=20;          % max hypotheses at each frame
p0.trms=0.03;                      % window length for rms measurement
p0.dtrms=0.02;                     % window spacing for rms measurement
p0.preemph=-7000;                  % s-plane position of preemphasis zero
p0.nfullag=7;                      % number of full lags to try (must be odd)
p=paramsetch(p0,q);

% redefine commonly used parameters

candtr=p.candtr;          % minimum peak in NCCF                      [0.3]
vtranc=p.vtranc;          % fixed voice-state transition cost         [0.005]
vtrac=p.vtrac;            % delta amplitude modulated transition cost [0.5]
vtrsc=p.vtrsc;            % delta spectrum modulated transition cost  [0.5]
vobias=p.vobias;          % bias to encourage voiced hypotheses       [0.0]
doublec=p.doublec;        % cost of exact doubling or halving         [0.35]
ncands=p.ncands;          % max hypotheses at each frame              [20]
nfullag=p.nfullag;        % number of full lags to try (must be odd)  [7]

% derived parameters (mostly dependent on sample rate fs)

krms=round(p.trms*fs);            % window length for rms measurement
kdrms=round(p.dtrms*fs);          % window spacing for rms measurement
rmswin=hanning(krms).^2;
kdsmp=round(0.25*fs/p.f0max);
hlpw=round(p.tlpw*fs/2);          % force window to be an odd length
blp=sinc((-hlpw:hlpw)/kdsmp).*hamming(2*hlpw+1).';
fsd=fs/kdsmp;
kframed=round(fsd*p.tframe);      % downsampled frame length
kframe=kframed*kdsmp;           % frame increment at full rate
rmsix=(1:krms)+floor((kdrms-kframe)/2); % rms index according to Talkin; better=(1:krms)+floor((kdrms-krms+1)/2)
minlag=ceil(fsd/p.f0max);
maxlag=round(fsd/p.f0min);        % use round() only because that is what Talkin does
kcorwd=round(fsd*p.tcorw);        % downsampled correlation window
kcorw=kcorwd*kdsmp;             % full rate correlation window
spoff=max(hlpw-floor(kdsmp/2),1+kdrms-rmsix(1)-kframe);  % offset for first speech frame at full rate
sfoff=spoff-hlpw+floor(kdsmp/2); % offset for downsampling filter
sfi=1:kcorwd;                   % initial decimated correlation window index array
sfhi=1:kcorw;                   % initial correlation window index array
sfj=1:kcorwd+maxlag;
sfmi=repmat((minlag:maxlag)',1,kcorwd)+repmat(sfi,maxlag-minlag+1,1);
lagoff=(minlag-1)*kdsmp;        % lag offset when converting to high sample rate
beta=p.lagwt*p.f0min/fs;            % bias towards low lags
log2=log(2);
lpcord=2+round(fs/1000);        % lpc order for itakura distance
hnfullag=floor(nfullag/2);
jumprat=exp((doublec+log2)/2);  % lag ratio at which octave jump cost is lowest
ssq=s.^2;
csssq=cumsum(ssq);
sqrt(min(csssq(kcorw+1:end)-csssq(1:end-kcorw))/kcorw);
afact=max([p.absnoise^2,max(ssq)*p.signoise^2,min(csssq(kcorw+1:end)-csssq(1:end-kcorw))*(p.relnoise/kcorw)^2])^2*kcorw^2;

% downsample signal to approx 2 kHz to speed up autocorrelation calculation
% kdsmp is the downsample factor

sf=filter(blp/sum(blp),1,s(sfoff+1:end));
sp=filter([1 exp(p.preemph/fs)],1,s); % preemphasised speech for LPC calculation
sf(1:length(blp)-1)=[];         % remove startup transient
sf=sf(1:kdsmp:end);             % downsample to =~2kHz
nsf=length(sf);                 % length of downsampled speech
ns=length(s);                   % length of full rate speech

% Calculate the frame limit to ensure we don't run off the end of the speech or decimated speech:
%   (a) For decimated autocorrelation when calculating sff():  (nframe-1)*kframed+kcorwd+maxlag <= nsf
%   (b) For full rate autocorrelation when calculating sfh():  max(fho)+kcorw+maxlag*kdsamp+hnfllag <= ns
%   (c) For rms ratio window when calculating rr            :  max(fho)+rmsix(end) <= ns
% where max(fho) = (nframe-1)*kframe + spoff

nframe=floor(1+min((nsf-kcorwd-maxlag)/kframed,(ns-spoff-max(kcorw-maxlag*kdsmp-hnfullag,rmsix(end)))/kframe));

% now search for autocorrelation peaks in the downsampled signal

cost=zeros(nframe,ncands);      % cumulative cost
prev=zeros(nframe,ncands);      % traceback pointer
mcands=zeros(nframe,1);         % number of actual candidates excluding voiceless
lagval=repmat(NaN,nframe,ncands-1);    % lag of each voiced candidate
tv=zeros(nframe,3);             % diagnostics: 1=voiceless cost, 2=min voiced cost, 3:cumulative voiceless-min voiced
if doback
    costms=cell(nframe,1);
end

% Main processing loop for each 10 ms frame

for iframe=1:nframe       % loop for each frame (~10 ms)
    
    % Find peaks in the normalized autocorrelation of subsampled (2Khz) speech
    % only keep peaks that are > 30% of highest peak
    
    fho=(iframe-1)*kframe+spoff;
    sff=sf((iframe-1)*kframed+sfj);
    sffdc=mean(sff(sfi));       % mean of initial correlation window length
    sff=sff-sffdc;              % subtract off the mean
    nccfd=normxcor(sff(1:kcorwd),sff(minlag+1:end));
    [ipkd,vpkd]=v_findpeaks(nccfd,'q');
    
    % Debugging: execute the line below to plot the autocorrelation peaks.
    % v_findpeaks(nccfd,'q'); xlabel(sprintf('Lag = (x+%d)*%g ms',minlag-1,1000*kdsmp/fs)); ylabel('Normalized Cross Correlation'); title (sprintf('Frame %d/%d',iframe,nframe));
    
    vipkd=[vpkd ipkd];
    vipkd(vpkd<max(vpkd)*candtr,:)=[];          % eliminate peaks that are small
    if size(vipkd,1)
        if size(vipkd,1)>ncands-1
            vipkd=sortrows(vipkd);
            vipkd(1:size(vipkd,1)-ncands+1,:)=[];   % eliminate lowest to leave only ncands-1
        end
        lagcan=round(vipkd(:,2)*kdsmp+lagoff);        % convert the lag candidate values to the full sample rate
        nlcan=length(lagcan);
    else
        nlcan=0;
    end
    
    % If there are any candidate lag values (nlcan>0) then refine their accuracy at the full sample rate
    
    if nlcan
        laglist=reshape(repmat(lagcan(:)',nfullag,1)+repmat((-hnfullag:hnfullag)',1,nlcan),nfullag*nlcan,1);
        sfh=s(fho+(1:kcorw+max(lagcan)+hnfullag));
        sfhdc=mean(sfh(sfhi));
        sfh=sfh-sfhdc;
        e0=sum(sfh(sfhi).^2);                     % energy of initial correlation window (only needed to store in tv(:,6)
        lagl2=repmat(lagcan(:)',nfullag+kcorw-1,1)+repmat((1-hnfullag:hnfullag+kcorw)',1,nlcan);
        nccf=normxcor(sfh(1:kcorw),sfh(lagl2),afact);
        
        [maxcc,maxcci]=max(nccf,[],1);
        vipk=[maxcc(:) lagcan(:)+maxcci(:)-hnfullag-1];
        vipk=vipk(:,[1 2 2]);
        maxccj=maxcci(:)'+nfullag*(0:nlcan-1);    % vector index into nccf array
        msk=mod(maxcci,nfullag-1)~=1 & 2*nccf(maxccj)-nccf(mod(maxccj-2,nfullag*nlcan)+1)-nccf(mod(maxccj,nfullag*nlcan)+1)>0;  % don't do quadratic interpolation for the end ones
        if any(msk)
            maxccj=maxccj(msk);
            vipk(msk,3)=vipk(msk,3)+(nccf(maxccj+1)-nccf(maxccj-1))'./(2*(2*nccf(maxccj)-nccf(maxccj-1)-nccf(maxccj+1)))';
        end
        vipk(maxcc<max(maxcc)*candtr,:)=[];          % eliminate peaks that are small
        if size(vipk,1)>ncands-1
            vipk=sortrows(vipk);
            vipk(1:size(vipk,1)-ncands+1,:)=[];   % eliminate lowest to leave only ncands-1
        end
        
        % vipk(:,1) has NCCF value, vipk(:,2) has integer peak position, vipk(:,3) has refined peak position
        
        mc=size(vipk,1);
    else
        mc=0;
    end
    
    % We now have mc lag candidates at the full sample rate
    
    mc1=mc+1;               % total number of candidates including "unvoiced" possibility
    mcands(iframe)=mc;      % save number of lag candidates (needed for pitch consistency cost calculation)
    if mc
        lagval(iframe,1:mc)=vipk(:,3)';
        cost(iframe,1)=vobias+max(vipk(:,1));   % voiceless cost
        cost(iframe,2:mc1)=1-vipk(:,1)'.*(1-beta*vipk(:,3)');   % local voiced costs
        tv(iframe,2)=min(cost(iframe,2:mc1));
    else
        cost(iframe,1)=vobias;          % if no lag candidates (mc=0), then the voiceless case is the only possibility
    end
    tv(iframe,1)=cost(iframe,1);
    if iframe>1                         % if it is not the first frame, then calculate pitch consistency and v/uv transition costs
        mcp=mcands(iframe-1);
        costm=zeros(mcp+1,mc1);         % cost matrix: rows and cols correspond to candidates in previous and current frames (incl voiceless)
        
        % if both frames have at least one lag candidate, then calculate a pitch consistency cost
        
        if mc*mcp
            lrat=abs(log(repmat(lagval(iframe,1:mc),mcp,1)./repmat(lagval(iframe-1,1:mcp)',1,mc)));
            costm(2:end,2:end)=p.freqwt*min(lrat,doublec+abs(lrat-log2));  % allow pitch doubling/halving
        end
        
        % if either frame has a lag candidate, then calculate the cost of voiced/voiceless transition and vice versa
        
        if mc+mcp
            rr=sqrt((rmswin'*s(fho+rmsix).^2)/(rmswin'*s(fho+rmsix-kdrms).^2)); % amplitude "gradient"
            ss=0.2/(distitar(lpcauto(sp(fho+rmsix),lpcord),lpcauto(sp(fho+rmsix-kdrms),lpcord),'e')-0.8);   % Spectral stationarity: note: Talkin uses Hanning instead of Hamming windows for LPC
            costm(1,2:end)= vtranc+vtrsc*ss+vtrac/rr;   % voiceless -> voiced cost
            costm(2:end,1)= vtranc+vtrsc*ss+vtrac*rr;
            tv(iframe,4:5)=[costm(1,mc1) costm(mcp+1,1)];
        end
        costm=costm+repmat(cost(iframe-1,1:mcp+1)',1,mc1);  % add in cumulative costs
        [costi,previ]=min(costm,[],1);
        cost(iframe,1:mc1)=cost(iframe,1:mc1)+costi;
        prev(iframe,1:mc1)=previ;
    else                            % first ever frame
        costm=zeros(1,mc1); % create a cost matrix in case doing a backward recursion
    end
    if mc
        tv(iframe,3)=cost(iframe,1)-min(cost(iframe,2:mc1));
        tv(iframe,6)=5*log10(e0*e0/afact);
    end
    if doback
        costms{iframe}=costm; % need to add repmatted cost into this
    end
end

% now do traceback

best=zeros(nframe,1);
[cbest,best(nframe)]=min(cost(nframe,1:mcands(nframe)+1));
for i=nframe:-1:2
    best(i-1)=prev(i,best(i));
end
vix=find(best>1);
fx=repmat(NaN,nframe,1);                        % unvoiced frames will be NaN
fx(vix)=fs*lagval(vix+nframe*(best(vix)-2)).^(-1); % leave as NaN if unvoiced
tt=zeros(nframe,3);
tt(:,1)=(1:nframe)'*kframe+spoff;       % find frame times
tt(:,2)=tt(:,1)+kframe-1;
jratm=(jumprat+1/jumprat)/2;
tt(2:end,3)=abs(fx(2:end)./fx(1:end-1)-jratm)>jumprat-jratm;    % new spurt if frequency ratio is outside (1/jumprat,jumprat)
tt(1,3)=1;           % first frame always starts a spurt
tt(1+find(isnan(fx(1:end-1))),3)=1; % NaN always forces a new spurt

% plot results if there are no output arguments of if the 'g' mode option is specified

if ~nargout || any(mode=='g')
    tf=spoff+(0:nframe-1)'*kframe;      % one sample before start of each frame
    blag=repmat(NaN,nframe,1);                        % unvoiced frames will be NaN
    blag(vix)=lagval(vix+nframe*(best(vix)-2)); % leave as NaN if unvoiced
    ts=(1:ns)/fs;                       % time scale for speech samples
    tsa=[1:tf(1) tf(end)+kframe+1:ns];  % indexes for unprocessed speech [-1 term is an error methinks]
    sup=repmat(NaN,ns,1);               % unprocessed speech - plot in black
    sup(tsa)=s(tsa);
    sv=reshape(s(tf(1)+1:tf(end)+kframe),kframe,nframe);               % processed speech
    su=sv;
    su(:,best>1)=NaN;                   % delete all voiced samples
    sv(:,best==1)=NaN;                  % delete all unvoiced samples
    tsuv=(tf(1)+1:tf(end)+kframe)/fs;
    su=su(:);
    sv=sv(:);
    ax=zeros(2,1);
    ax(1)=subplot(211);
    plot(ts,sup,'-k',tsuv,su,'r-',tsuv,sv,'b-');
    title('Speech');
    ax(2)=subplot(212);
    plot((tf+(kframe+1)/2)/fs,lagval*1000/fs,'xr',(tf+(kframe+1)/2)/fs,blag*1000/fs,'-b')
    xlabel('Time (s)');
    ylabel('Period (ms)');
    title('Lag Candidates');
    linkaxes(ax,'x');
end
if ~any(mode=='u')
    tt(isnan(fx),:)=[];    % remove NaN spurts
    fx(isnan(fx),:)=[];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v=normxcor(x,y,d)
% Calculate the normalized cross correlation of column vectors x and y
% we can calculate this in two ways but fft is much faster even for nx small
% We must have nx<=ny and the output length is ny-nx+1
% note that this routine does not do mean subtraction even though this is normally a good idea
% if y is a matrix, we correlate with each column
% d is a constant added onto the normalization factor
% v(j)=x'*yj/sqrt(d + x'*x * yj'*yj) where yj=y(j:j+nx-1) for j=1:ny-nx+1

if nargin<3
    d=0;
end
nx=length(x);
[ny,my]=size(y);
nv=1+ny-nx;
if nx>ny
    error('second argument is shorter than the first');
end

nf=pow2(nextpow2(ny));
w=irfft(repmat(conj(rfft(x,nf,1)),1,my).*rfft(y,nf,1));
s=zeros(ny+1,my);
s(2:end,:)=cumsum(y.^2,1);
v=w(1:nv,:)./sqrt(d+(x'*x).*(s(nx+1:end,:)-s(1:end-nx,:)));