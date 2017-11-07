function [d,g,rr,ss]=sigalign(s,r,maxd,m,fs)
%SIGALIGN align a clean reference with a noisy signal [d,g,rr,ss]=(s,r,maxd,m,fs)
% Inputs:
%            m  mode
%                 u = unity gain
%                 g = find optimal gain [default]
%                 a = A-weight the signals
%                 b = weight signals by BS-468
%                 s = find delay to maximize the correlation coefficient between r and s [default]
%                 S = find delay to maximize the energy of the component of r in s
%                 p = plot result
%            s  test signal
%            r  reference signal
%         maxd  [+-max] or [min max] delay allowed in samples or fractions of length(r)
%               default is maximum that ensures at least 50% of r or s in the overlap
%           fs  sample frequency (only used for filtering and plotting)
%
% Outputs:
%            d = optimum delay to apply to r
%            g = optimal gain to apply to r
%           rr = g*r(* -d)  [zero padded to match s if ss output is not given]
%           ss = s truncated if necessary to martch to the length of rr


%      Copyright (C) Mike Brookes 2011
%      Version: $Id: sigalign.m 713 2011-10-16 14:45:43Z dmb $
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

% Bugs/Suggestions
% 1. add option to calculate a DC offset
% 2. optionally find optimal fractional time shift
% 3. split long signals into chunks to reduce memory requirements

ns=length(s);
nr=length(r);
if numel(s)~=ns || numel(r)~=nr
    error('Inputs cannot be matrices');
end
s=s(:);
r=r(:);
if nargin<3
    maxd=[];
end
switch numel(maxd)
    case 0
        if nr<ns
            lmm=[-0.25*nr ns-0.75*nr];
        else
            lmm=[-0.25*ns nr-0.75*ns];
        end
    case 1
        lmm=[-maxd maxd];
    otherwise
        lmm=maxd(1:2);
end
lmm=round(lmm.*(1+(nr-1)*(abs(lmm)<1)));  % convert fractions of nr to samples
lmin=lmm(1);
lmax=lmm(2);
lags=lmax-lmin+1;
if lags<=0
    error('Invalid lag limits');
end
if nargin<4 || ~numel(m)
    m='gs';
end
if nargin<5 || ~numel(fs)
    fs=[];
else
    if any(m=='a')
        [b,a]=stdspectrum(2,'z',fs);
        s=filter(b,a,s);
        r=filter(b,a,r);
    elseif any(m=='b')
        [b,a]=stdspectrum(8,'z',fs);
        s=filter(b,a,s);
        r=filter(b,a,r);
    end
end

% now do cross correlation

rxi=max(1,1-lmin);   % first reference sample needed
rxj=min(nr,ns-lmax); % last reference sample needed
nrx=rxj-rxi+1;          % length of reference segment
if nrx<1
    error('Reference signal too short');
end
fl=2^nextpow2(lmax-lmin+nrx);
sxi=max(1,rxi+lmin); % first signal sample needed
sxj=min(ns,rxj+lmax); % last signal sample needed
rs=irfft(rfft([s(sxi:sxj); zeros(fl-sxj+sxi-1,1)]).*conj(rfft([r(rxi:rxj); zeros(fl-rxj+rxi-1,1)])));
rsu=rs(1:lags);
ssq=cumsum(s(sxi:sxj).^2);
ssqd=[ssq(nrx); ssq(nrx+1:lmax-lmin+nrx)-ssq(1:lmax-lmin)];
if any (m=='S') % maximize energy of common component
    [cmx,icx]=max(abs(rsu)); % maximize cross correlation
else
    [cmx,icx]=max(rsu.^2./ssqd); % maximize correlation coefficient
end
d=icx-1+lmin;
ia=max(1,d+1); % first sample of s in common region
ja=min(ns,d+nr); % last sample of s in common region
ija=ia:ja;
ijad=ija-d;
rr=r(ijad);
ss=s(ija);
if any (m=='u')
    g=1;
else
g=sum(rr.*ss)/sum(rr.^2);   % gain to apply to r
end
rr=rr*g;
if ~nargout || any(m=='p')
    xco=sum(rr.*ss)/sqrt(sum(rr.^2)*sum(ss.^2));
    snr=sum(rr.^2)/sum((rr-ss).^2);
    if numel(fs)==1
        tun='s';
    else
        tun='samples';
        fs=1;
    end
    subplot(311);
    plot(ija/fs,rr);
    pm='+-';
    title(sprintf('Ref delay = %.2g %s, %cGain = %.2g dB, Xcorr = %.2g, SNR = %.2g dB',d/fs,tun,pm(1+(g<0)),20*log10(g),xco,10*log10(snr)));
    ylabel('Reference');
    set(gca,'XLim',ija([1 end])/fs);
    axh(2)=gca;
    subplot(312);
    plot(ija/fs,ss);
    ylabel('Signal');
    set(gca,'XLim',ija([1 end])/fs);
    axh(1)=gca;
    subplot(313);
    plot(ija/fs,ss-rr);
    ylabel('Residual');
    xlabel(sprintf('Time (%s)',tun));
    set(gca,'XLim',ija([1 end])/fs);
    axh(3)=gca;
    linkaxes(axh(1:3),'x');
end
if nargout==3
    rr=[zeros(ia-1,1); rr; zeros(ns-ja,1)]; % force to be the size of s
end



