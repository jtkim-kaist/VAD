function [b,a,fx,bx,gd]=gammabank(n,fs,w,fc,bw,ph,k)
%GAMMABANK gammatone filter bank [b,a,fx,bx,gd]=(n,fs,w,fc,bw,ph,k)
%
% Usage:
%          (1) [b,a,fx,bx,gd]=gammabank(0.35,fs,'',[100 6000]);
%              y = filterbank(b,a,s,gd);
%
%              Will create an erb-spaced filterbank between 100 Hz and 6kHz
%              with a filter spacing of 0.35 erb and a default bandwidth
%              of 1.019 erb. Omitting the "y =" from the second line will plot
%              a spectrogram.
% Inputs:
%       n   number of filters in filterbank or the filter spacing in
%           bark/erb/k-mel. Set n=0 if fc lists centre frequencies explicitly.
%		fs  sample rate in Hz
%		fc  centre frequencies [default = [100 6000] ]
%		bw  bandwidths [default = 1.019*erb(fc) ]
%       ph  phase of gammatone impulse repsonse [default = 0]
%       k   filter order [default = 4]
%		w   any sensible combination of the following:
%             'e' = erb scale for filter spacing, frequencies ('F' option)
%             and bandwidths ('W' option but can be overridden by 'EMBLH')
%                   'm','b','l','h' = mel, bark, log10 and Hz scale
%             'E' = erb scale for bandwidths ('W' option)
%                   'M','B','L','H' = mel, bark, log10 and Hz scale
%
%             'n' = n input gives number of filters [default if n>=1]
%             'N' = n input gives filter spacing  [default if n<1]
%
%             'f' = fc is in Hz [default]
%             'F' = fc is in mel/erb-rate/bark/log10
%             'w' = bw inputs are in Hz [default]
%             'W' = bw inputs are in multiples of df/dx where x=bark/erb/mel etc
%
%             'k' = force a filter at 1kHz
%             ['d' = choose ph() so that all filters have zero DC gain]
%             ['a' = use all-pole gammtone funtion: see [1]]
%             ['s' = use Slaney gammatone approximation: see [2]]
%             ['z' = use one-zero gammatone function: see [1]]
%
%             'g' = plot filter responses [default if no output arguments present]
%             'G' = plot frequency responses on a log axis
%
% Outputs:
%       b/a    filter coefficients: one filter per row
%       fx,bx  centre frequencies and bandwidths in Hz
%       gd     group delay at the centre frequencies (in samples)
%

% The impulse response of filter i is proportional to:
%       h(n)=((n/fs).^(k-1))*cos(2*pi*fx(i)*n/fs+ph(i))*exp(-2*pi*bx(i)*n/fs)
% where n=0,1,2,...
% Note that the DC gain is only equal to zero for one particular value of ph(i)
% The filters are normalized to have unity gain at the centre frequencies
%
% References
%  [1]	R. F. Lyon, A. G. Katsiamis, and E. M. Drakakis.
%       History and future of auditory filter models.
%       In Proc Intl Symp Circuits and Systems, pages 3809–3812, 2010.
%       doi: 10.1109/ISCAS.2010.5537724.
%  [2]	M. Slaney.
%       An efficient implementation of the patterson-holdsworth auditory filter bank.
%       Technical report, Apple Computer, Perception Group, Tech. Rep, 1993.

%      Copyright (C) Mike Brookes 2009-2010
%      Version: $Id: gammabank.m 7847 2016-04-29 06:55:47Z dmb $
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

if nargin<7
    k=[];
    if nargin<6
        ph=[];
        if nargin<5
            bw=[];
            if nargin<4
                fc=[];
                if nargin<3
                    w='';
                end
            end
        end
    end
end
fx=fc(:);
bx=bw(:);
if ~numel(k)
    k=4;
end
if ~numel(fx)
    fx=[100; 6000]; % default
end
wr='e';   % default frequency warping is erb
for i=1:length(w)
    if any(w(i)=='bmlef');
        wr=w(i);
    end
end
if any(w=='k')
    fk=1000;
    switch wr              % convert 1kHz to spacing units
        case 'b'
            fk=frq2bark(fk);
        case 'l'
            fk=log10(fk);
        case 'e'
            fk=frq2erb(fk);
    end
else
    fk=0;
end
if any(w=='W')
    wb=wr;
else
    wb='h';     % default bandwidth units are Hz
end
for i=1:length(w)
    if any(w(i)=='BMLEF');
        wb=w(i)+'a'-'A';        % convert to lower case
    end
end
if ~numel(bx)
    bx=1.019;
    wb='e';
end
if any(w=='F')          % convert centre frequencies to Hz
    switch wr
        case 'b'
            fx=bark2frq(fx);
        case 'm'
            fx=mel2frq(fx);
        case 'l'
            fx=10.^(fx);
        case 'e'
            fx=erb2frq(fx);
    end
end

% now sort out the centre frequencies

if n>0                      % n>0: filter end points specified
    bx=bx(1);               % only use the first bx element
    if n==1             % only one filter requested
        fx=fx(1);           % just use the first frequency
    else
        switch wr               % convert end frequencies to spacing units
            case 'b'
                fx=frq2bark(fx);
            case 'm'
                fx=frq2mel(fx);
            case 'l'
                fx=log10(fx);
            case 'e'
                fx=frq2erb(fx);
        end
        if n<1 || any(w=='N')       % n = filter spacing
            if fk               % force filter to 1 kHz
                f0=fk-n*floor((fk-fx(1))/n);
            else                % centre filters in range
                f0=(fx(2)+fx(1)-n*floor((fx(2)-fx(1))/n))/2;
            end
            fx=(f0:n:fx(2))';

        else                        % n = number of filters specified
            % Multiple filters - evenly spaced
            fx=linspace(fx(1),fx(2),n)';     % centre frequencies in spacing units
            if fk              % force a filter at 1kHz
                ik=1+ceil((fk-fx(1))*(n-1)/(fx(n)-fx(1))); % index of centre freq immediately above 1 kHz
                if ik>n || ik>1 && ((fk-fx(1))*(fx(n)-fx(ik-1))>(fx(n)-fk)*(fx(ik)-fx(1)))
                    fx=fx(1)+(fx-fx(1))*(fk-fx(1))/(fx(ik)-fx(1));
                else
                    fx=fx(n)+(fx-fx(n))*(fx(n)-fk)/(fx(n)-fx(ik-1));
                end
            end
        end
        switch wr % convert back to Hz
            case 'b'
                fx=bark2frq(fx);
            case 'm'
                fx=mel2frq(fx);
            case 'l'
                fx=10.^(fx);
            case 'e'
                fx=erb2frq(fx);
        end
    end

end
% now sort out the bandwidths
nf=numel(fx);
if numel(bx)==1
    bx=bx(ones(nf,1));      % replicate if necessary
end
switch wb               % convert bandwidth to Hz
    case 'b'
        [dum,bwf]=frq2bark(fx);
    case 'm'
        [dum,bwf]=frq2mel(fx);
    case 'l'
        bwf=fx*log(10);
    case 'e'
        [dum,bwf]=frq2erb(fx);
    case 'h'
        bwf=ones(nf,1);
end
bx=bx.*bwf;
if ~numel(ph)
    ph=0;
end
if numel(ph)==1
    ph=ph(ones(nf,1));      % replicate if necessary
else
    ph=ph(:);
end
%
% t=(0:ceil(10*fs/(2*pi*bnd)))/fs;  % five time constants
% gt=t.^(n-1).*cos(2*pi*cfr*t+phi).*exp(-2*pi*bnd*t);
% gt=gt/sqrt(mean(gt.^2)); % normalize
% figure(1);
% plot(t,gt);
% title('Desired Impulse response');
% xlabel(['Time (' xticksi 's)']);
%
ww=exp((1i*fx-bx)*2*pi/fs);             % pole position in top half of z-plane
a=round([1 cumprod((-k:-1)./(1:k))]);   % create binomial coefficients
b=conv(a,(0:k-1).^(k-1));               % convolve with powers of k-1
b=exp(1i*ph)*b(1:k);                    % correct for starting phase shift
wwp=repmat(ww,1,k+1).^repmat(0:k,nf,1); % powers of pole position
denc=repmat(a,nf,1).*wwp;               % replicate pole k times
numc=b.*wwp(:,1:k);
b=zeros(nf,2*k);                        % space for numerators
a=zeros(nf,2*k+1);                      % space for denominators
gd=zeros(nf,1);                         % space for group delay
ww=exp(2i*fx*pi/fs);                    % exp(j*centre-freq)
for i=1:nf
    b(i,:)=real(conv(numc(i,:),conj(denc(i,:))));
    a(i,:)=real(conv(denc(i,:),conj(denc(i,:)))); % denominator has k repeated poles at ww and ww'
    u=polyval(b(i,:),ww(i));
    v=polyval(a(i,:),ww(i));
    ud=polyval(b(i,:).*(0:2*k-1),ww(i));
    vd=polyval(a(i,:).*(0:2*k),ww(i));
    b(i,:)=b(i,:)*abs(v/u);
    gd(i)=real((v*ud-u*vd)/(u*v));     % group delay at centre freq in samples
end

% now plot graph

if ~nargout || any(w=='g') || any(w=='G')
    ng=200;      %number of points to plot
    if any(w=='G')
        fax=logspace(log10(fx(1)/4),log10(fs/2),ng);
    else
        fax=linspace(0,fs/2,ng);
    end
    ww=exp(2i*pi*fax/fs);
    gg=zeros(nf,ng);
    for i=1:nf
        gg(i,:)=10*log10(abs(polyval(b(i,:),ww)./polyval(a(i,:),ww)));
    end
    if any(w=='G')
        semilogx(fax,gg','-b');
        set(gca,'xlim',[fax(1) fax(end)]);
    else
        plot(fax,gg','-b');
    end
    xlabel(['Frequency (' xticksi 'Hz)']);
    set(gca,'ylim',[-50 1]);
    title(sprintf('Order-%d Gammatone Filterbank (N=%d, Opt=%s)',k,nf,w));
end
