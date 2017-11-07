function [x,zo]=estnoiseg(yf,tz,pp)
%ESTNOISEG - estimate MMSE noise spectrum [x,zo]=(yf,tz,pp)
%
% Usage:    ninc=round(0.016*fs);   % frame increment [fs=sample frequency]
%           ovf=2;                  % overlap factor
%           f=rfft(enframe(s,hanning(ovf*ninc,'periodic'),ninc),ovf*ninc,2);
%           f=f.*conj(f);           % convert to power spectrum
%           x=estnoiseg(f,ninc/fs); % estimate the noise power spectrum
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
%
% The algorithm parameters are defined in reference [1] from which equation
% numbers are given in parentheses. They are as follows:
%
%        pp.tax      % smoothing time constant for noise power estimate [0.0717 seconds](8)
%        pp.tap      % smoothing time constant for smoothed speech prob [0.152 seconds](23)
%        pp.psthr    % threshold for smoothed speech probability [0.99] (24)
%        pp.pnsaf    % noise probability safety value [0.01] (24)
%        pp.pspri    % prior speech probability [0.5] (18)
%        pp.asnr     % active SNR in dB [15] (18)
%        pp.psini    % initial speech probability [0.5] (23)
%        pp.tavini   % assumed speech absent time at start [0.064 seconds]
%
% If convenient, you can call estnoiseg in chunks of arbitrary size. Thus the following are equivalent:
%
%                   (a) dp=estnoiseg(yp(1:300),tinc);
%
%                   (b) [dp(1:100,:),z]=estnoiseg(yp(1:100,:),tinc);
%                       [dp(101:200,:),z]=estnoiseg(yp(101:200,:),z);
%                       [dp(201:300,:),z]=estnoiseg(yp(201:300,:),z);


% This is intended to be a precise implementation of [1] for a frame rate of 62.5 Hz.
% Time constants are adjusted for other frame rates.
%
% Refs:
%    [1] Gerkmann, T. & Hendriks, R. C.
%        Unbiased MMSE-Based Noise Power Estimation With Low Complexity and Low Tracking Delay
%        IEEE Trans Audio, Speech, Language Processing, 2012, 20, 1383-1393

%	   Copyright (C) Mike Brookes 2012
%      Version: $Id: estnoiseg.m 3387 2013-08-23 12:32:47Z dmb $
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
if isempty(yf) && isstruct(tz)             % no real data
    zo=tz;                  % just keep the same state
else
    if isstruct(tz)         % take parameters from a previous call
        nrcum=tz.nrcum;     % cumulative number of frames
        xt=tz.xt;           % smoothed power spectrum
        pslp=tz.pslp;       % correction factor (9)
        tinc=tz.tinc;       % frame increment
        qq=tz.qq;           % parameter structure
    else
        tinc = tz;          % second argument is frame increment
        nrcum=0;            % no frames so far
        % default algorithm constants
        qq.tax=0.0717;      % noise output smoothing time constant = -tinc/log(0.8) (8)
        qq.tap=0.152;       % speech prob smoothing time constant = -tinc/log(0.9) (23)
        qq.psthr=0.99;      % threshold for smoothed speech probability [0.99] (24)
        qq.pnsaf=0.01;      % noise probability safety value [0.01] (24)
        qq.pspri=0.5;       % prior speech probability [0.5] (18)
        qq.asnr=15;         % active SNR in dB [15] (18)
        qq.psini=0.5;       % initial speech probability [0.5] (23)
        qq.tavini=0.064;        % assumed speech absent time at start [64 ms]

        if nargin>=3 && ~isempty(pp)  % update fields from pp input
            qqn=fieldnames(qq);
            for i=1:length(qqn)
                if isfield(pp,qqn{i})
                    qq.(qqn{i})=pp.(qqn{i});
                end
            end
        end
        pslp=repmat(qq.psini,1,nrf); % initialize smoothed speech presence prob
        xt=[];                       % initialize just in case the first call has no data
    end

    % unpack parameters needed within the loop

    psthr=qq.psthr;     % threshold for smoothed speech probability [0.99] (24)
    pnsaf=qq.pnsaf;     % noise probability safety value [0.01] (24)

    % derived algorithm constants

    ax=exp(-tinc/qq.tax); % noise output smoothing factor = 0.8 (8)
    axc=1-ax;
    ap=exp(-tinc/qq.tap); % noise output smoothing factor = 0.9 (23)
    apc=1-ap;
    xih1=10^(qq.asnr/10); % speech-present SNR
    xih1r=1/(1+xih1)-1;
    pfac=(1/qq.pspri-1)*(1+xih1); % p(noise)/p(speech) (18)

    if nrcum==0 && nr>0       % initialize values for first frame
        xt=qq.psini*mean(yf(1:max(1,min(nr,round(1+qq.tavini/tinc))),:),1);  % initial noise estimate
    end

    % loop for each frame
    for t=1:nr
        yft=yf(t,:);        % noisy speech power spectrum
        ph1y=(1+pfac*exp(xih1r*yft./xt)).^(-1); % a-posteriori speech presence prob (18)
        pslp=ap*pslp+apc*ph1y; % smoothed speech presence prob (23)
        ph1y=min(ph1y,1-pnsaf*(pslp>psthr)); % limit ph1y (24)
        xtr=(1-ph1y).*yft+ph1y.*xt; % estimated raw noise spectrum (22)
        xt=ax*xt+axc*xtr;  % smooth the noise estimate (8)
        x(t,:)=xt;  % save the noise estimate
    end
    if nargout>1    % we need to store the state for next time
        zo.nrcum=nrcum+nr;      % number of frames so far
        zo.xt=xt;          % smoothed power spectrum
        zo.pslp=pslp;               % correction factor (9)
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

