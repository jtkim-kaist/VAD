function [x,zo,xs]=estnoisem(yf,tz,pp)
%ESTNOISEM - estimate noise spectrum using minimum statistics
%
% Usage:    ninc=round(0.016*fs);   % frame increment [fs=sample frequency]
%           ovf=2;                  % overlap factor
%           f=rfft(enframe(s,hanning(ovf*ninc,'periodic'),ninc),ovf*ninc,2);
%           f=f.*conj(f);           % convert to power spectrum
%           x=estnoisem(f,ninc/fs); % estimate the noise power spectrum
%
% Inputs:
%   yf      input power spectra (one row per frame)
%   tz      frame increment in seconds
%           Alternatively, the input state from a previous call (see below)
%   pp      algorithm parameters [optional]
%
% Outputs:
%   x       estimated noise power spectra (one row per frame)
%   zo      output state
%   xs      estimated std error of x (one row per frame)
%           xs seems often to be an underestimate by a factor of 2 or 3
%
% The algorithm parameters are defined in reference [1] from which equation
% numbers are given in parentheses. They are as follows:
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
% Example use:      y=enframe(s,w,ni);                  % divide speech signal s(n) into
%                                                       % overlapping frames using window w(n)
%                   yf=rfft(y,nf,2);                    % take fourier transform
%                   dp=estnoisem(yf.*conj(yf),tinc);    % estimate the noise
%
% If convenient, you can call estnoisem in chunks of arbitrary size. Thus the following are equivalent:
%
%                   (a) dp=estnoisem(yp(1:300),tinc);
%
%                   (b) [dp(1:100),z]=estnoisem(yp(1:100),tinc);
%                       [dp(101:200),z]=estnoisem(yp(101:200),z);
%                       [dp(201:300),z]=estnoisem(yp(201:300),z);


% This is intended to be a precise implementation of [1] with Table III
% replaced by the updated table 5 from [2]. The only deliberate algorithm
% change is the introduction of a minimum value for 1/Qeq in equation (23).
% This change only affects the first few frames and improves the
% convergence of the algorithm. A minor improveemnt was reported in [3] but
% this has not yet been included.
%
% Refs:
%    [1] Rainer Martin.
%        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
%        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.
%    [2] Rainer Martin.
%        Bias compensation methods for minimum statistics noise power spectral density estimation
%        Signal Processing, 2006, 86, 1215-1229
%    [3] Dirk Mauler and Rainer Martin
%        Noise power spectral density estimation on highly correlated data
%        Proc IWAENC, 2006

%	   Copyright (C) Mike Brookes 2008
%      Version: $Id: estnoisem.m 1718 2012-03-31 16:40:41Z dmb $
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

[nr,nrf]=size(yf);          % number of frames and freq bins
x=zeros(nr,nrf);            % initialize output arrays
xs=zeros(nr,nrf);           % will hold std error in the future
if isempty(yf) && isstruct(tz)             % no real data
    zo=tz;              % just keep the same state
else
    if isstruct(tz)       % take parameters from a previous call
        nrcum=tz.nrcum;
        p=tz.p;          % smoothed power spectrum
        ac=tz.ac;               % correction factor (9)
        sn2=tz.sn2;              % estimated noise power
        pb=tz.pb;               % smoothed noisy speech power (20)
        pb2=tz.pb2;
        pminu=tz.pminu;
        actmin=tz.actmin;   % Running minimum estimate
        actminsub=tz.actminsub;           % sub-window minimum estimate
        subwc=tz.subwc;                   % force a buffer switch on first loop
        actbuf=tz.actbuf;  % buffer to store subwindow minima
        ibuf=tz.ibuf;
        lminflag=tz.lminflag;      % flag to remember local minimum
        tinc=tz.tinc;     % frame increment
        qq=tz.qq;         % parameter structure
    else
        tinc = tz;          % second argument is frame increment
        nrcum=0;            % no frames so far
        % default algorithm constants

        qq.taca=0.0449;    % smoothing time constant for alpha_c = -tinc/log(0.7) in equ (11)
        qq.tamax=0.392;    % max smoothing time constant in (3) = -tinc/log(0.96)
        qq.taminh=0.0133;    % min smoothing time constant (upper limit) in (3) = -tinc/log(0.3)
        qq.tpfall=0.064;   % time constant for P to fall (12)
        qq.tbmax=0.0717;   % max smoothing time constant in (20) = -tinc/log(0.8)
        qq.qeqmin=2;       % minimum value of Qeq (23)
        qq.qeqmax=14;      % max value of Qeq per frame
        qq.av=2.12;             % fudge factor for bc calculation (23 + 13 lines)
        qq.td=1.536;       % time to take minimum over
        qq.nu=8;           % number of subwindows
        qq.qith=[0.03 0.05 0.06 Inf]; % noise slope thresholds in dB/s
        qq.nsmdb=[47 31.4 15.7 4.1];

        if nargin>=3 && ~isempty(pp)
            qqn=fieldnames(qq);
            for i=1:length(qqn)
                if isfield(pp,qqn{i})
                    qq.(qqn{i})=pp.(qqn{i});
                end
            end
        end
    end

    % unpack parameter structure

    taca=qq.taca;    % smoothing time constant for alpha_c = -tinc/log(0.7) in equ (11)
    tamax=qq.tamax;    % max smoothing time constant in (3) = -tinc/log(0.96)
    taminh=qq.taminh;    % min smoothing time constant (upper limit) in (3) = -tinc/log(0.3)
    tpfall=qq.tpfall;   % time constant for P to fall (12)
    tbmax=qq.tbmax;   % max smoothing time constant in (20) = -tinc/log(0.8)
    qeqmin=qq.qeqmin;       % minimum value of Qeq (23)
    qeqmax=qq.qeqmax;      % max value of Qeq per frame
    av=qq.av;             % fudge factor for bc calculation (23 + 13 lines)
    td=qq.td;       % time to take minimum over
    nu=qq.nu;           % number of subwindows
    qith=qq.qith; % noise slope thresholds in dB/s
    nsmdb=qq.nsmdb;   % maximum permitted +ve noise slope in dB/s

    % derived algorithm constants

    aca=exp(-tinc/taca); % smoothing constant for alpha_c in equ (11) = 0.7
    acmax=aca;          % min value of alpha_c = 0.7 in equ (11) also = 0.7
    amax=exp(-tinc/tamax); % max smoothing constant in (3) = 0.96
    aminh=exp(-tinc/taminh); % min smoothing constant (upper limit) in (3) = 0.3
    bmax=exp(-tinc/tbmax); % max smoothing constant in (20) = 0.8
    snrexp = -tinc/tpfall;
    nv=round(td/(tinc*nu));    % length of each subwindow in frames
    if nv<4            % algorithm doesn't work for miniscule frames
        nv=4;
        nu=max(round(td/(tinc*nv)),1);
    end
    nd=nu*nv;           % length of total window in frames
    [md,hd]=mhvals(nd); % calculate the constants M(D) and H(D) from Table III
    [mv,hv]=mhvals(nv); % calculate the constants M(D) and H(D) from Table III
    nsms=10.^(nsmdb*nv*tinc/10);  % [8 4 2 1.2] in paper
    qeqimax=1/qeqmin;  % maximum value of Qeq inverse (23)
    qeqimin=1/qeqmax; % minumum value of Qeq per frame inverse

    if isempty(yf)      % provide dummy initialization
        ac=1;               % correction factor (9)
        subwc=nv;                   % force a buffer switch on first loop
        ibuf=0;
        p=x;          % smoothed power spectrum
        sn2=p;              % estimated noise power
        pb=p;               % smoothed noisy speech power (20)
        pb2=pb.^2;
        pminu=p;
        actmin=repmat(Inf,1,nrf);   % Running minimum estimate
        actminsub=actmin;           % sub-window minimum estimate
        actbuf=repmat(Inf,nu,nrf);  % buffer to store subwindow minima
        lminflag=zeros(1,nrf);      % flag to remember local minimum
    else

        if ~nrcum       % initialize values for first frame
            p=yf(1,:);          % smoothed power spectrum
            ac=1;               % correction factor (9)
            sn2=p;              % estimated noise power
            pb=p;               % smoothed noisy speech power (20)
            pb2=pb.^2;
            pminu=p;
            actmin=repmat(Inf,1,nrf);   % Running minimum estimate
            actminsub=actmin;           % sub-window minimum estimate
            subwc=nv;                   % force a buffer switch on first loop
            actbuf=repmat(Inf,nu,nrf);  % buffer to store subwindow minima
            ibuf=0;
            lminflag=zeros(1,nrf);      % flag to remember local minimum
        end

        % loop for each frame

        for t=1:nr              % we use t instead of lambda in the paper
            yft=yf(t,:);        % noise speech power spectrum
            acb=(1+(sum(p)./sum(yft)-1).^2).^(-1);  % alpha_c-bar(t)  (9)
            ac=aca*ac+(1-aca)*max(acb,acmax);       % alpha_c(t)  (10)
            ah=amax*ac.*(1+(p./sn2-1).^2).^(-1);    % alpha_hat: smoothing factor per frequency (11)
            snr=sum(p)/sum(sn2);
            ah=max(ah,min(aminh,snr^snrexp));       % lower limit for alpha_hat (12)

            p=ah.*p+(1-ah).*yft;            % smoothed noisy speech power (3)
            b=min(ah.^2,bmax);              % smoothing constant for estimating periodogram variance (22 + 2 lines)
            pb=b.*pb + (1-b).*p;            % smoothed periodogram (20)
            pb2=b.*pb2 + (1-b).*p.^2;     	% smoothed periodogram squared (21)

            qeqi=max(min((pb2-pb.^2)./(2*sn2.^2),qeqimax),qeqimin/(t+nrcum));   % Qeq inverse (23)
            qiav=sum(qeqi)/nrf;             % Average over all frequencies (23+12 lines) (ignore non-duplication of DC and nyquist terms)
            bc=1+av*sqrt(qiav);             % bias correction factor (23+11 lines)
            bmind=1+2*(nd-1)*(1-md)./(qeqi.^(-1)-2*md);      % we use the simplified form (17) instead of (15)
            bminv=1+2*(nv-1)*(1-mv)./(qeqi.^(-1)-2*mv);      % same expression but for sub windows
            kmod=bc*p.*bmind<actmin;        % Frequency mask for new minimum
            if any(kmod)
                actmin(kmod)=bc*p(kmod).*bmind(kmod);
                actminsub(kmod)=bc*p(kmod).*bminv(kmod);
            end
            if subwc>1 && subwc<nv              % middle of buffer - allow a local minimum
                lminflag=lminflag | kmod;    	% potential local minimum frequency bins
                pminu=min(actminsub,pminu);
                sn2=pminu;
            else
                if subwc>=nv                    % end of buffer - do a buffer switch
                    ibuf=1+rem(ibuf,nu);     	% increment actbuf storage pointer
                    actbuf(ibuf,:)=actmin;    	% save sub-window minimum
                    pminu=min(actbuf,[],1);
                    i=find(qiav<qith);
                    nsm=nsms(i(1));          	% noise slope max
                    lmin=lminflag & ~kmod & actminsub<nsm*pminu & actminsub>pminu;
                    if any(lmin)
                        pminu(lmin)=actminsub(lmin);
                        actbuf(:,lmin)=repmat(pminu(lmin),nu,1);
                    end
                    lminflag(:)=0;
                    actmin(:)=Inf;
                    subwc=0;
                end
            end
            subwc=subwc+1;
            x(t,:)=sn2;
            qisq=sqrt(qeqi);
            % empirical formula for standard error based on Fig 15 of [2]
            xs(t,:)=sn2.*sqrt(0.266*(nd+100*qisq).*qisq/(1+0.005*nd+6/nd)./(0.5*qeqi.^(-1)+nd-1));
        end
    end
    if nargout>1    % we need to store the state for next time
        zo.nrcum=nrcum+nr;      % number of frames so far
        zo.p=p;          % smoothed power spectrum
        zo.ac=ac;               % correction factor (9)
        zo.sn2=sn2;              % estimated noise power
        zo.pb=pb;               % smoothed noisy speech power (20)
        zo.pb2=pb2;
        zo.pminu=pminu;
        zo.actmin=actmin;   % Running minimum estimate
        zo.actminsub=actminsub;           % sub-window minimum estimate
        zo.subwc=subwc;                   % force a buffer switch on first loop
        zo.actbuf=actbuf;  % buffer to store subwindow minima
        zo.ibuf=ibuf;
        zo.lminflag=lminflag;      % flag to remember local minimum
        zo.tinc=tinc;     % must be the last one
        zo.qq=qq;
    end
    if ~nargout
        clf;
        subplot(212);
        plot((1:nr)*tinc,10*log10([sum(yf,2) sum(x,2)]))
        ylabel('Frame Energy (dB)');
        xlabel(sprintf('Time (s)   [%d ms frame incr]',round(tinc*1000)));
        axisenlarge([-1 -1.05]);
        legend('input','noise','Location','Best');
        subplot(211);
        plot(1:nrf,10*log10([sum(yf,1)'/nr sum(x,1)'/nr]))
        ylabel('Power (dB)');
        xlabel('Frequency bin');
        axisenlarge([-1 -1.05]);
        legend('input','noise','Location','Best');
    end
end

function [m,h,d]=mhvals(d)
% Values are taken from Table 5 in [2]
%[2] R. Martin,"Bias compensation methods for minimum statistics noise power
%               spectral density estimation", Signal Processing Vol 86, pp1215-1229, 2006.

% approx: plot(d.^(-0.5),[m 1-d.^(-0.5)],'x-'), plot(d.^0.5,h,'x-')
persistent dmh
if isempty(dmh)
    dmh=[
        1   0       0;
        2   0.26    0.15;
        5   0.48    0.48;
        8   0.58    0.78;
        10  0.61    0.98;
        15  0.668   1.55;
        20  0.705   2;
        30  0.762   2.3;
        40  0.8     2.52;
        60  0.841   3.1;
        80  0.865   3.38;
        120 0.89    4.15;
        140 0.9     4.35;
        160 0.91    4.25;
        180 0.92    3.9;
        220 0.93    4.1;
        260 0.935   4.7;
        300 0.94    5];
end

if nargin>=1
    i=find(d<=dmh(:,1));
    if isempty(i)
        i=size(dmh,1);
        j=i;
    else
        i=i(1);
        j=i-1;
    end
    if d==dmh(i,1)
        m=dmh(i,2);
        h=dmh(i,3);
    else
        qj=sqrt(dmh(i-1,1));    % interpolate using sqrt(d)
        qi=sqrt(dmh(i,1));
        q=sqrt(d);
        h=dmh(i,3)+(q-qi)*(dmh(j,3)-dmh(i,3))/(qj-qi);
        m=dmh(i,2)+(qi*qj/q-qj)*(dmh(j,2)-dmh(i,2))/(qi-qj);
    end
else
    d=dmh(:,1);
    m=dmh(:,2);
    h=dmh(:,3);
end