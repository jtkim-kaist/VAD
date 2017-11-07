function [lev,xx] = activlevg(sp,fs,mode)
%ACTIVLEVG Measure active speech level robustly [LEV,AF,FSO]=(sp,FS,MODE)
%
%Inputs: sp     is the speech signal
%        FS     is the sample frequency in Hz (see also FSO below)
%        MODE   is a combination of the following:
%               r - raw omit input filters (default is 200 Hz to 5.5 kHz)
%               0 - no high pass filter (i.e. include DC)
%               4 - high pass filter at 40 Hz (but allows mains hum to pass)
%               1 - use cheybyshev 1 filter
%               2 - use chebyshev 2 filter (default)
%               e - use elliptic filter
%               h - omit low pass filter at 5.5 kHz
%               d - give outputs in dB rather than power
%               n - output a normalized speech signal as the first argument
%               N - output a normalized filtered speech signal as the first argument
%               l - give additional power level estimates (see below for details)
%               a - include A-weighting filter
%               i - include ITU-R-BS.468/ITU-T-J.16 weighting filter
%
%Outputs:
%    If the "n" option is specified, a speech signal normalized to 0dB will be given as
%    the first output followed by the other outputs.
%        LEV    gives the speech level in units of power (or dB if mode='d')
%               if mode='l' is specified, LEV is a row vector containing:
%                       [ active-level mean-power mean-noise-power P.56-level harmonic-power-level]

% This is an implementation of the algorithm described in [1].
%
% Refs:
%    [1] S. Gonzalez and M. Brookes.
%        Speech active level estimation in noisy conditions.
%        In Proc. IEEE Intl Conf. Acoustics, Speech and Signal Processing,
%        pp 6684–6688, Vancouver, May 2013. doi: 10.1109/ICASSP.2013.6638955.

%      Copyright (C) Sira Gonzalez, Mike Brookes 2008-2012
%      Version: $Id: activlevg.m 4968 2014-08-05 18:22:14Z dmb $
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
persistent ro snr_ro offset c25zp c15zp e5zp
if isempty(ro)
    ro = [0,0.124893,0.160062,0.212256,0.279864,0.351281,0.437356,...
        0.550314,0.680075,0.795741,0.892972,0.939087,1];  %  mapping from SNR to rho
    snr_ro = -2:0.5:4;
    offset = 0.8472; % correction factor for harmonic power in dB
    c25zp=[0.37843443673309i 0.23388534441447i; -0.20640255179496+0.73942185906851i -0.54036889596392+0.45698784092898i];
    c25zp=[[0; -0.66793268833792] c25zp conj(c25zp)];
    %       [c1z,c1p,c1k] = cheby1(5,0.25,1,'high','s');
    c15zp=[-0.659002835294875+1.195798636925079i -0.123261821596263+0.947463030958881i];
    c15zp=[zeros(1,5); -2.288586431066945 c15zp conj(c15zp)];
    %      [ez,ep,ek] = ellip(5,0.25,50,1,'high','s')
    e5zp=[0.406667680649209i 0.613849362744881i; -0.538736390607201+1.130245082677107i -0.092723126159100+0.958193646330194i];
    e5zp=[[0; -1.964538608244084]  e5zp conj(e5zp)];
    %    w=linspace(0.2,2,100);
    %    figure(1); plot(w,20*log10(abs(freqs(real(poly(c15zp(1,:))),real(poly(c15zp(2,:))),w)))); title('Chebyshev 1');
    %    figure(2); plot(w,20*log10(abs(freqs(real(poly(c25zp(1,:))),real(poly(c25zp(2,:))),w)))); title('Chebyshev 2');
    %    figure(3); plot(w,20*log10(abs(freqs(real(poly(e5zp(1,:))),real(poly(e5zp(2,:))),w)))); title('Elliptic');
end
if nargin<3
    mode=' ';
end
fso.ffs=fs;                       	% sample frequency
if any(mode=='r')                   % included for backward compatibility
    mode=['0h' mode];               % abolish both filters
elseif fs<14000
    mode=['h' mode];               % abolish lowpass filter at low sample rates
end
% s-plane zeros and poles of high pass 5'th order filter -0.25dB at w=1 and -50dB stopband
if any(mode=='1')
    szp=c15zp;            % Chebyshev 1
elseif any(mode=='e')
    szp=e5zp;             % Elliptic
else
    szp=c25zp;            % Chebyshev 2
end
% calculate high pass filter at 40 or 200 Hz
if all(mode~='0')
    if any(mode=='4')
        fl=40;               % 40 Hz cutoff
    else
        fl=200;              % 200 Hz cutoff
    end
    zl=2./(1-szp*tan(fl*pi/fs))-1;      % 40 or 200 Hz LF limit
    al=real(poly(zl(2,:)));          	% high pass filter
    bl=real(poly(zl(1,:)));
    sw=(-1).^(0:5)';                    %z^(-n) for z=-1
    bl=bl*(al*sw)/(bl*sw);       	% scale to give HF gain of 1
end
% calculate low pass filter at 5500 Hz
if all(mode~='h')
    zh=2./(szp/tan(5500*pi/fs)-1)+1;
    ah=real(poly(zh(2,:)));
    bh=real(poly(zh(1,:)));
    bh=bh*sum(ah)/sum(bh);
end
if any(mode=='a')
    [bw aw]=stdspectrum(2,'z',fs);
elseif any(mode=='i')
    [bw aw]=stdspectrum(8,'z',fs);
end
% apply the input filters to the speech
if all(mode~='0')
    sq=filter(bl,al,sp(:));     % highpass filter
else
    sq=sp(:);
end
if all(mode~='h')
    sq=filter(bh,ah,sq(:));     % lowpass filter
end
if any(mode=='a') || any(mode=='i')
    sq=filter(bw,aw,sq(:));     % weighting filter
end
p56 = activlev(sq,fs,'0hl');            % get active level from P.56 method (no extra filters)
[fx,dum,pv]=fxpefac(sp,fs);              % Estimate f0 and voiced probability from unfiltered speech
[tz,f,S]=spgrambw(sq,fs,20,[0 5 fs/2],[],0.01);   % spectrogram with 20 Hz bandwidth and 5 Hz/10 ms resolution

nh = 15;         % Number of harmonics to evaluate
tv = find(pv>0.5);  % Find voiced frames
nv=length(tv);      % number of voiced frames
ran = -60:60;       % Frequency range of the mexican hat (+/- 60Hz)
nm = length(ran);   % width of mexican hat
o = 15;             % semi-width of central lobe
mexic = (1-(ran/o).^2).*exp(-(ran.^2)/(2*o.^2)); % Mexican hat wavelet
%% calculate harmonic energy in voiced frames
Et=zeros(nv,1);
for jv =1:nv                                    % loop for each voiced frame
    f_har = fx(tv(jv)).*(1:nh).';               % Calculate frequencies of harmonics
    lev = repmat(f_har,1,nm)+repmat(ran,nh,1);  % frequencies to test: nh x nm
    data = interp1(f,S(tv(jv),:),lev,'spline'); % interpolate spectrogram onto required frequencies
    Eh = sum(data.*repmat(mexic,nh,1),2);       % estimate power in each harmonic
    Et(jv) = sum(Eh(Eh>0));                     % sum the non-negative harmonic powers
end
har = mean(Et);                                 % mean power of voiced frames (not including offset)
%% Noise estimate
dp=estnoiseg(S,tz(2)-tz(1));                        % estimate noise power spectrum
dpm=mean((f(2)-f(1))*sum(dp,2));                    % mean noise psd
snr_est=10*log10(har/dpm);  % estimate global SNR (should this be restricted to voiced frames?)
%% Combine methods

levp56 = 10*log10(p56(1));
levharm = 10*log10(har)+offset;
if snr_est<=-2
    lev = levharm;
elseif snr_est>4
    lev = levp56;
else
    ro1 = interp1(snr_ro,ro,snr_est);               % interpolate the value of rho
    lev = ro1*levp56 + (1-ro1)*levharm;
end
levdb=lev;
if any(mode=='l')
%                       [ active-level mean-power mean-noise-power P.56-level mean-harmonic-power]
    lev=[levdb 10*log10(p56(2)) 10*log10(dpm) levp56 levharm];
end
if all(mode~='d')
    lev=10.^(0.1*lev);
end
if any(mode=='n') || any(mode=='N')
    xx=lev;
    if any(mode=='n')
        sq=sp;
    end
    if dpm>0
        lev=sq*10^(-0.05*levdb);
    else
        lev=sq;
    end
end