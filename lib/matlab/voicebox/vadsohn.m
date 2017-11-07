function [vs,zo]=vadsohn(si,fsz,m,pp)
%VADSOHN implements a voice activity detector [VS,ZO]=(S,FSZ,M,P)
%
% Inputs:
%   si      input speech signal
%   fsz     sample frequency in Hz
%           Alternatively, the input state from a previous call (see below)
%   m       mode [default = 'a']:
%               'n'   output frame start/end in samples
%               't'   output frame start/end in seconds
%               'a'   output activity decision [default]
%               'b'   output activity likelihood ratio
%               'p'   plot graph [default if no outputs]
%   pp      algorithm parameters [optional]
%
% Outputs:
%   vs(:,n)     outputs in the order [t1,t2,a,b] as selected by m options
%               if 'n' or 't' is specified in m, vs has one row per frame and
%               t1 and t2 give the frame start end end times (in samples or seconds).
%               otherwise, vs has one row per sample.
%   zo          output state (used as input for a subsequent call).
%
% The algorithm operation is controlled by a small number of parameters:
%
%        pp.of          % overlap factor = (fft length)/(frame increment) [2]
%        pp.ti          % desired output frame increment (10 ms)
%        pp.tj          % internal frame increment (10 ms)
%        pp.ri          % set to 1 to round tj to the nearest power of 2 samples [0]
%        pp.ta          % Time const for smoothing SNR estimate [0.396 seconds]
%        pp.gx          % maximum posterior SNR as a power ratio [1000 = 30dB]
%        pp.gz          % minimum posterior SNR as a power ratio [0.0001 = -40dB]
%        pp.xn          % minimum prior SNR [0]
%        pp.pr          % Speech probability threshold [0.7]
%        pp.ts          % mean talkspurt length (100 ms)
%        pp.tn          % mean silence length (50 ms)
%        pp.ne          % noise estimation: 0=min statistics, 1=MMSE [0]
%
% In addition it is possible to specify parameters for the noise estimation algorithm
% which implements reference [3] from which equation numbers are given in parentheses.
% They are as follows:
%
%        pp.taca      % (11): smoothing time constant for alpha_c [0.0449 seconds]
%        pp.tamax     % (3): max smoothing time constant [0.392 seconds]
%        pp.taminh    % (3): min smoothing time constant (upper limit) [0.0133 seconds]
%        pp.tpfall    % (12): time constant for P to fall [0.064 seconds]
%        pp.tbmax     % (20): max smoothing time constant [0.0717 seconds]
%        pp.qeqmin    % (23): minimum value of Qeq [2]
%        pp.qeqmax    % max value of Qeq per frame [14]
%        pp.av        % (23)+13 lines: fudge factor for bc calculation  [2.12]
%        pp.td        % time to take minimum over [1.536 seconds]
%        pp.nu        % number of subwindows to use [3]
%        pp.qith      % Q-inverse thresholds to select maximum noise slope [0.03 0.05 0.06 Inf ]
%        pp.nsmdb     % corresponding noise slope thresholds in dB/second   [47 31.4 15.7 4.1]
%
%
% If convenient, you can call vadsohn in chunks of arbitrary size. Thus the following are equivalent:
%
%                   (a) y=vadsohn(s,fs);
%
%                   (b) [y1,z]=vadsohn(s(1:1000),fs);
%                       [y2,z]=vadsohn(s(1001:2000),z);
%                       y3=vadsohn(s(2001:end),z);
%                       y=[y1; y2; y3];
%
% Note that in all cases the number of output samples will be less than the number of input samples if
% the input length is not an exact number of frames. In addition, if two output arguments
% are specified, the last partial frame samples will be retained for overlap adding
% with the output from the next call to specsub().
%
% Refs:
%    [1] J. Sohn, N. S. Kim, and W. Sung.
%        A statistical model-based voice activity detection.
%        IEEE Signal Processing Lett., 6 (1): 1?, 1999. doi: 10.1109/97.736233.
%    [2] Ephraim, Y. & Malah, D.
%        Speech enhancement using a minimum-mean square error short-time spectral amplitude estimator
%        IEEE Trans Acoustics Speech and Signal Processing, 32(6):1109-1121, Dec 1984
%    [3] Rainer Martin.
%        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
%        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.

%      Copyright (C) Mike Brookes 2004
%      Version: $Id: vadsohn.m 3100 2013-06-13 16:05:56Z dmb $
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
if nargin<3 || isempty(m)
    m='';
end
if isstruct(fsz)
    fs=fsz.fs;
    qq=fsz.qq;
    qp=fsz.qp;
    ze=fsz.ze;
    s=zeros(length(fsz.si)+length(si(:)),1); % allocate space for speech
    s(1:length(fsz.si))=fsz.si;
    s(length(fsz.si)+1:end)=si(:);
else
    fs=fsz;     % sample frequency
    s=si(:);
    % default algorithm constants

    qq.of=2;        % overlap factor = (fft length)/(frame increment)
    qq.pr=0.9;      % Speech probability threshold
    qq.ts=0.1;      % mean talkspurt length (100 ms)
    qq.tn=0.05;     % mean silence length (50 ms)
    qq.ti=10e-3;    % desired output frame increment (10 ms)
    qq.tj=10e-3;    % internal frame increment (10 ms)
    qq.ri=0;        % round ni to the nearest power of 2
    qq.ta=1000;    % Time const for smoothing SNR estimate = -tinc/log(0.98) from [2]
    qq.gx=1000;     % maximum posterior SNR = 30dB
    qq.gz=1e-4;     % minimum posterior SNR = -40dB
    qq.xn=0;        % minimum prior SNR = -Inf dB
    qq.ne=1;          % noise estimation: 0=min statistics, 1=MMSE [0]
    if nargin>=4 && ~isempty(pp)
        qp=pp;      % save for estnoisem call
        qqn=fieldnames(qq);
        for i=1:length(qqn)
            if isfield(pp,qqn{i})
                qq.(qqn{i})=pp.(qqn{i});
            end
        end
    else
        qp=struct;  % make an empty structure
    end
end
% derived algorithm constants
qq.tj=min([qq.tj 0.5*qq.ts 0.5*qq.tn]);     % at least two frames per talk/silence spurt
nj=max(round(qq.ti/qq.tj),1);               % number of internal frames per output frame
if qq.ri
    ni=pow2(nextpow2(qq.ti*fs*sqrt(0.5)/nj)); % internal frame increment in samples
else
    ni=round(qq.ti*fs/nj);    % internal frame increment in samples
end
tinc=ni/fs;             % true frame increment time
a=exp(-tinc/qq.ta);     % SNR smoothing coefficient
gmax=qq.gx;             % max posterior SNR = 30 dB
kk=sqrt(2*pi);          % sqrt(8)*Gamma(1.5) - required constant
xn=qq.xn;            	% floor for prior SNR, xi
gz=qq.gz;               % floor for posterior SNR
a01=tinc/qq.tn;         % a01=P(silence->speech)
a00=1-a01;              % a00=P(silence->silence)
a10=tinc/qq.ts;         % a10=P(speech->silence)
a11=1-a10;              % a11=P(speech->speech)
b11=a11/a10;
b01=a01/a00;
b00=a01-a00*a11/a10;
b10=a11-a10*a01/a00;
prat=a10/a01;       % P(silence)/P(speech)
lprat=log(prat);
% calculate power spectrum in frames

no=round(qq.of);                                   % integer overlap factor
nf=ni*no;           % fft length
nd=floor(ni*(no-1)/2); % output delay relative to start of frame
w=hamming(nf+1)'; w(end)=[];  % for now always use hamming window
w=w/sqrt(sum(w(1:ni:nf).^2));      	% normalize to give overall gain of 1
ns=length(s);
y=enframe(s,w,ni);
yf=rfft(y,nf,2);
if ~size(yf,1)                                  % no data frames
    vs=[];
    nr=0;
    nb=0;
else
    yp=yf.*conj(yf);                % power spectrum of input speech
    [nr,nf2]=size(yp);              % number of frames
    nb=ni*nr;
    isz=isstruct(fsz);
    if isz
        if qq.ne>0
            [dp,ze]=estnoiseg(yp,ze);	% estimate the noise using MMSE
        else
            [dp,ze]=estnoisem(yp,ze);	% estimate the noise using minimum statistics
        end
        xu=fsz.xu;                  % previously saved unsmoothed SNR
        lggami=fsz.gg;
        nv=fsz.nv;
    else
        if qq.ne>0
            [dp,ze]=estnoiseg(yp,tinc,qp);	% estimate the noise using MMSE
        else
            [dp,ze]=estnoisem(yp,tinc,qp);	% estimate the noise using minimum statistics
        end
        xu=1;                               % dummy unsmoothed SNR from previous frame [2](53)++
        lggami=0;                            % initial prob ratio
        nv=0;                               % number of previous speech samples
    end
    gam=max(min(yp./dp,gmax),gz);       % gamma = posterior SNR [2](10)
    prsp=zeros(nr,1);                   % create space for prob ratio vector
    for i=1:nr                          % loop for each frame
        gami=gam(i,:);                  % gamma(i) = a posteriori SNR [2](10)
        xi=a*xu+(1-a)*max(gami-1,xn);   % xi = smoothed a priori SNR [2](53)
        xi1=1+xi;
        v=0.5*xi.*gami./xi1;            % note that this is 0.5*vk in [2]
        lamk=2*v-log(xi1);              % log likelihood ratio for each frequency bin [1](3)
        lami=sum(lamk(2:nf2))/(nf2-1);  % mean log LR over frequency omitting DC term [1](4)
        if lggami<0                     % avoid overflow in calculating [1](11)
            lggami=lprat+lami+log(b11+b00/(a00+a10*exp(lggami)));
        else
            lggami=lprat+lami+log(b01+b10/(a10+a00*exp(-lggami)));
        end
        if lggami<0
            gg=exp(lggami);
            prsp(i)=gg/(1+gg);
        else
            prsp(i)=1/(1+exp(-lggami));
        end
        gi=(0.277+2*v)./gami;           % accurate to 0.02 dB for v>0.5
        mv=v<0.5;
        if any(mv)                      % only calculate Bessel functions for v<0.5
            vmv=v(mv);
            gi(mv)=kk*sqrt(vmv).*((0.5+vmv).*besseli(0,vmv)+vmv.*besseli(1,vmv))./(gami(mv).*exp(vmv)); % [2](7)
        end
        xu=gami.*gi.^2;                 % unsmoothed prior SNR % [2](53)
    end
    ncol=any(m=='a')+any(m=='b');       % number of output data columns
    if ~ncol
        m(end+1)='a';   % force at least one output
        ncol=1;
    end
    nw=floor(nr/nj);                % number of output or plotted frames
    if any(m=='n') || any(m=='t')       % frame-based outputs
        vs=zeros(nw,2+ncol);            % space for outputs
        vs(:,1)=nd+1+nv+(0:nw-1)'*ni*nj;   % starting sample of each frame
        vs(:,2)=ni*nj-1+vs(:,1);           % ending sample of each frame
        if any(m=='t')
            vs(:,1:2)=vs(:,1:2)/fs;     % convert to seconds
        end
        if any(m=='a')
            vs(:,3)=any(reshape(prsp(1:nw*nj)>qq.pr,nj,[]),1).';
        end
        if any(m=='b')
            vs(:,end)=max(reshape(prsp(1:nw*nj),nj,[]),[],1).';
        end
    else
        na=nd*(1-isz);                  % include preamble if no input state supplied
        nc=(nargout<=1)*(ns-nd-nb);     % include tail if no output state desired
        vs=zeros(na+nb+nc,ncol);
        vs(na+(1:nb),ncol)=reshape(repmat(prsp',ni,1),nb,1);
        vs(1:na,ncol)=vs(na+1,ncol);        % replicate start
        vs(na+nb+1:end,ncol)=vs(na+nb,ncol); % replicate end
        if any(m=='a')
            vs(:,1)=vs(:,ncol)>qq.pr;
        end
    end
end
if nargout>1
    zo.si=s(nd+nb+1:end);      % save the tail end of the input speech signal
    zo.fs=fs;                       % save sample frequency
    zo.qq=qq;                       % save local parameters
    zo.qp=qp;                       % save estnoisem parameters
    zo.ze=ze;                       % save state of noise estimation
    zo.xu=xu;                       % unsmoothed prior SNR  [2](53)
    zo.gg=lggami;                    % posterior prob ratio: P(speech|s)/P(silence|s) [1](11)
    zo.nv=nv+nd+nb;                      % number of previous speech samples
end
if (~nargout || any(m=='p')) && nr>0
    ax=subplot(212);
    plot((nv+nd+[0 nr*ni])/fs,[qq.pr qq.pr],'r:',(nv+nd+ni*nj/2+(0:nw-1)*ni*nj)/fs,max(reshape(prsp(1:nw*nj),nj,[]),[],1).','b-' );
    set(gca,'ylim',[-0.05 1.05]);
    xlabel('Time (s)');
    ylabel('Pr(speech)');
    ax(2)=subplot(211);
    plot((nv+(1:ns))/fs,s);
    ylabel('Speech');
    title('Sohn VAD');
    linkaxes(ax,'x');
end