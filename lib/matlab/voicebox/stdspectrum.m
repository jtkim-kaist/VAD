function [b,a,si,sn]=stdspectrum(s,m,f,n,zi,bs,as)
%STDSPECTRUM Generate standard acoustic/speech spectra in s- or z-domain [B,A,SI,SN]=(S,M,F,N,ZI,BS,AS)
%
% Usage: (1) [b,a]=stdspectrum(2,'z',fs)  % create A-weighting z-domain filter for sample frequency fs
%        (2) [b,a]=stdspectrum(2,'zMLT',fs) % as above but also plot both s- and z-domain responses with log freq axis
%        (3) x=stdspectrum(11,'t',fs,n) % generate n samples of speech-shaped noise
%        (4) [b,a]=stdspectrum(0,'zEMLT',fs,[],[],bs,as) % convert s-domain filter bs(s)/as(s) to z-domain and plot approximation error
%        (5) for i=1:10; figure(i); stdspectrum(i,'z',3e4); end; tilefigs; % plot all the spectra for fs=30kHz
%
%Inputs:  s  Spectrum type (either text or number - see below) or 0 to use bs/as
%         m  mode: char 1 specifies output type (no combinations allowed),
%                    f - frequency response (complex)
%                    m - magnitude response
%                    p - power spectrum
%                    l - power per decade
%                    d - decibel power spectrum
%                    e - decibel power per decade
%                    t - time waveform
%                    s - s-domain filter: b(s)/a(s) [default]
%                    z - z-domain filter: b(z)/a(z)
%                    i - sampled impulse response
%                plotting options
%                    M - plot magnitude spectrum in dB
%                    E - plot magnitude spectrum error in dB
%                    F - magnitude spectrum error in dB with a -40 dB floor relative to peak
%                    Q - plot phase spectrum
%                    T - also plot target spectra
%                    A - plot zeros/poles
%                    U - plot using a uniform frequency axis
%                    L - plot using a log frequency axis
%                    W - waveform
%                    S - spectrogram
%         f  sample frequency (modes 'z','i','t') or set of frequencies in Hz (modes 'f','m','p','d')
%         n  number of output samples (mode 'i','t')
%         zi initial state of filter from a previous call (mode 't')
%         bs numerator s-domain polynomial (or cell array containing polynomial factors) if s=0
%         as denominator s-domain polynomial (or cell array containing polynomial factors if s=0
%
% Outputs:  b       (1) numerator of the output spectrum (modes 's' or 'z')
%                   (2) output waveform (mode 't')
%                   (3) outptut spectrum (modes 'f', 'm', 'p' or 'd')
%           a       (1) denonminator of the output spectrum (modes 's' or 'z')
%                   (2) final state of the filter - use as the zi input of a future call (mode 't')
%           si     spectrum type number (0 to 10)
%           sn     spectrum name
%
% Spectrum types (specify either as a number or case-insensitive text abbreviation):
%   0  external     : BS and AS arguments specify an s-domain filter
%   1  White        : white noise
%   2  A-Weight     : the formula for this is given in [3] and is based on
%                     the equal-loudness curves of [9]
%   3  B-Weight     : this needs to be confirmed with ANSI S1.4-1981 standard or IEC 60651
%   4  C-Weight     : the formula for this is given in [3]
%   7  SII-intinv   : The inverse spectrum of the ear's internal masking noise; this is taken
%                     from table 1 of [1]. It is inverted so that it is a bandpass rather than
%                     bandstop characteristic.
%   8  BS-468       : The weighting proposed for audio frequency noise measurement in [5] and [6].
%   9  USASI        : Noise simulating long-term programme material spectrum from [7],[8].
%                     The level is such that the power is 0dB over an infinite bandwidth
%  10  POTS         : the D spectrum from [11].
%  11  LTASS-P50    : the long-term average speech spectrum taken from Table 1 of [4].
%                     Converted from mouth reference point @ 0.025m to dB SPL @ 1m on-axis.
%  13  LTASS-1994   : the long-term average speech spectrum that is taken from Table 2 in [2]
%
% Obsolete fits included for backward compatibility only:
%
%   5  X1-LTASS-P50  : (use 11 instead) the long-term average speech spectrum taken from Table 1 of [4].
%   6  X1-LTASS-1994 : (use 13 instead) the long-term average speech spectrum that is taken from Table 2 in [2]
%  12  X2-LTASS-1994 : (use 13 instead) the long-term average speech spectrum that is taken from Table 2 in [2]

% References:
% [1]	Methods for the calculation of the speech intelligibility index.
%       ANSI Standard S3.5-1997 (R2007), American National Standards Institute, 1997.
% [2]	D. Byrne, H. Dillon, K. Tran, S. Arlinger, K. Wilbraham, R. Cox, B. Hayerman,
%       R. Hetu, J. Kei, C. Lui, J. Kiessling, M. N. Kotby, N. H. A. Nasser,
%       W. A. H. E. Kholy, Y. Nakanishi, H. Oyer, R. Powell, D. Stephens, R. Meredith,
%       T. Sirimanna, G. Tavartkiladze, G. I. Frolenkov, S. Westerman, and C. Ludvigsen.
%       An international comparison of long-term average speech spectra.
%       JASA, 96 (4): 2108–2120, Oct. 1994.
% [3]	CENELEC. Electroacoustics - sound level meters. Technical Report EN EN 61672-1:2003, 2003.
%       (also ANSI S1.42-2001)
% [4]	ITU-T. Artificial voices. Standard P.50, Sept. 1999.
% [5]   ITU-T. Measurement of weighted noise in sound-programme circuits.
%       Recommendation J.16, 1988.
% [6]   ITU-R. Measurement of audio-requency noise voltage level in sound broadcasting.
%       Recommendation BS.468., 1986.
% [7]   NRSC AM Reemphasis, Deemphasize, and Broadcast Audio Transmission Bandwidth Specifications,
%       EIA-549 Standard, Electronics Industries Association , July 1988.
% [8]   NRSC AM Reemphasis, Deemphasize, and Broadcast Audio Transmission Bandwidth Specifications,
%       NRSC-1-A Standard, Sept 2007, Online: http://www.nrscstandards.org/SG/NRSC-1-A.pdf
% [9]   H. Fletcher and W. A. Munson. Loudness, its definition, measurement and calculation.
%       J. Acoust Soc Amer, 5: 82–108, Oct. 1933.
% [10]  American National Standard Specification for Sound Level Meters.
%       ANSI S1.4-1983 (R2006)/ANSI S1.4a-1985 (R2006), American National Standards Institute
% [11]	IEEE standard equipment requirements and measurement techniques for analog transmission
%       parameters for telecommunications. Standard IEEE Std 743-1995, Dec. 1995.

% Other candidates: (a) Z-weighting, (b) ISO226, (c) P.48 spectra
%
% Other standards:
%    IEEE743 has several weighting filters defined
%    ITU-T 0.41 Psophometer for use on telephone-type circuits
%    Bell System Technical Reference 41009 (C-message)
%    ISO 8041:2005 (E): Human Response to Vibration – Measuring
%    Instrumentation
%    IEC 1260:1995, class 1 (also IEC 61260/ANSI S1.11-2004) Octave band and fractional octave band filters
%    IEC 651: Specification for Sound Level Meters
%    IRS P.48: sending and receiving characteristics defined by isolated points
%    mIRS P.830 modified IRS also defined by isolated points (see Annex D) available in G.191
%    G.191 software tools library contains IRS and mIRS implementations in FIR and IIR

%      Copyright (C) Mike Brookes 2008
%      Version: $Id: stdspectrum.m 8211 2016-07-20 20:59:16Z dmb $
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

% Bugs/Suggestions:
% * Could generate impulse response by appending a LP filter (elliptic) to
%   the s-domain transfer function and sampling the impulse response
% * better calculation of impulse response length based on its total power
% * ensure that the number of z-domain zeros at z=1 is correct.

persistent spty nspty spz ient fixz baz bazps
% spty contains the name of the spectrum
% spz contains a list of the poles and zeros
%    spz(1) = gain, spz(2) = number of zeros (excluding implicit conjugates), spz(3...) zeros followed by poles
%    only one of a complex conjugate pair is given. To print out b/a in the correct format use:
%        ra=roots(a);ra=ra(imag(ra)>=0);rb=roots(b);rb=rb(imag(rb)>=0);fprintf('\n[%0.14g %d',b(1)/a(1),length(rb));
%        for i=1:length(rb), fprintf(' %s',sprintcpx(rb(i),'0.14gi')); end
%        for i=1:length(ra), fprintf(' %s',sprintcpx(ra(i),'0.14gi')); end, fprintf('];\n');
if isempty(spz)
    spty={'White';'A-Weight';'B-Weight';'C-Weight';'X1-LTASS-P50';'X1-LTASS-1994';'SII-IntInv';'BS-468';'USASI';'POTS';'LTASS-P50';'X2-LTASS-1994';'LTASS-1994'};
    nspty=size(spty,1);
    spz={[1 0];
        [7390100803.6603 4 0 0 0 0 -129.42731565506 -129.42731565506 -676.40154023295 -4636.125126885 -76618.526016858 -76618.526016858];
        [5986190155.0156 3 0 0 0 -129.42731565506 -129.42731565506 -995.88487118796 -76618.526016858 -76618.526016858];
        [5912384617.784 2 0 0 -129.42731565506 -129.42731565506 -76618.526016858 -76618.526016858];
        [1.1294790345421e+015 3 0 0 -34437.856184098 -721.94747118664+754.97798119504i -1721.704402273 -5234.2950286868 -10953.570773216+42789.342252749i];
        [19720493.192959 5 0 0 0 -22550.637578954 -11319.635610404+70239.177107659i -253.31327846696+672.10855509952i -1299.1885437118+2301.2064056419i -10646.952627978+68290.702816027i -147307.51763333];
        [6.1311266394354e+018 2 0 0 -381.08293630892 -5920.974797779 -4701.76218192+24369.279310049i 10597.854874768+39258.915154617i];
        [2.1034520039796e+024 1 0 -25903.701047817 -23615.535213635+36379.908937329i -62675.170058468 -18743.746690721+62460.156452506i];
        [72.648989380657 1 0 -2*pi*[100 320]];
        [7.8820088171767e+016 4 0 0 0 0 -452.681+1924.28i -2334+1702.73i -11264.2+8213.32i -4665.8+19828.7i];
        [2.67e+013 3 0 0 -34437.856184098 -721.94747118664+754.97798119504i -1721.704402273 -5234.2950286868 -10953.570773216+42789.342252749i];
        [1600000000 5 0 0 0 -22550.637578954 -11319.635610404+70239.177107659i -253.31327846696+672.10855509952i -1299.1885437118+2301.2064056419i -10646.952627978+68290.702816027i -147307.51763333  -628.3185307180 ];
        [0.40294607181247 6 0 0 0 0 -103821.41495527 -6138.6617784378 -1387.7780129857+2694.482976041i -212.70001558505+701.9845877834i -489.69089908576 -241.48780111882];
        };
    bazps=cell(nspty,4);
    for i=1:nspty
        spzi=spz{i};
        nsz=spzi(2);
        sz=spzi(3:3+nsz-1);
        sz=[sz conj(sz(imag(sz)~=0))];
        sp=spzi(3+nsz:end);
        sp=[sp conj(sp(imag(sp)~=0))];
        bazps{i,1}=spzi(1)*poly(sz);
        bazps{i,2}=poly(sp);
        bazps{i,3}=sz;
        bazps{i,4}=sp;
    end
    nz=15;  % size of cache
    ient=0; % cache entry number
    fixz=repmat(-1,nz,2);
    baz=cell(nz,2);

end
if nargin<2 || ~numel(m)
    m=' ';
end
m1=m(1);        % output format
if ~any(m1=='fmpldetszi')
    m1=char('s'+~nargout*('d'-'s')); % 's' normally, 'd' if no outputs
end
if nargin<3
    f=8192;  % default frequency
end

% determine the spectrum type

if ~numel(s) || s(1)==0
    si=0;
    sn='';
    sb=1;
    sz=[];  % list of s-domain zeros
    if iscell(bs)
        for i=1:numel(bs)
            sb=conv(sb,bs{i});
            sz=[sz roots(bs{i}).'];
        end
    elseif numel(bs)>0
        sb=bs;
        sz=roots(bs).';
    end
    sa=1;
    sp=[];  % list of s-domain poles
    if iscell(as)
        for i=1:numel(as)
            sa=conv(sa,as{i});
            sp=[sp roots(as{i}).'];
        end
    elseif numel(as)>0
        sa=as;
        sp=roots(as).';
    end
else
    if ischar(s)
        si=find(strcmpi(s,spty));
        if isempty(si)
            error('undefined spectrum type: %s',s);
        end
    else
        si=s;
    end
    if si>nspty
        error('undefined spectrum type: %d',si);
    end
    sn=spty{si};              % name of spectrum
    % get s-domain function
    sb=bazps{si,1};  % sb/sa is transfer function
    sa=bazps{si,2};
    sz=bazps{si,3};  % sz and sp are lists of zeros and poles
    sp=bazps{si,4};
end
if (nargin<3 || ~numel(f)) && any(m1=='fmpd') % calcualte the frequency range
    apz=abs([sp sz]);
    apz(apz==0)=[]; % ignore zero frequency poles/zeros
    if ~numel(apz)
        apz=[100 5000];
    elseif length(apz)==1
        apz=[apz/10 apz*10];
    end
    f=logspace(log10(min(apz)*0.5/pi)-0.5,log10(max(apz)*0.5/pi)+0.5);
end
if any(m1=='fmpdle')
    h=polyval(sb,2i*pi*f)./polyval(sa,2i*pi*f);
end
fs=f; % save sampling frequency
if any(m1=='izt') % we need a z-domain filter
    if si==1 % treat white noise specially
        bz=1;
        az=1;
    else
        jent=find(fixz(:,2)==f*nspty+si,1); % see if it is in the cache already
        if numel(jent)
            ient=ient+1; % cache entry number
            fixz(jent,1)=ient;  % update index to show it is recently used
            bz=baz{jent,1};
            az=baz{jent,2};
        else
            % we use an iterative method to find the best digital filter
            % we initialize the phases with either bilinear or impulse invariance
            % only using the impulse invariance if it is good (max error < 10dB)
            % we then iterate invfreqz using the s-domain magnitudes and
            % the phases of the best fit so far.
            % we use log-spaced frequencies at low frequencies and linear at high
            % we then search for various numbers of poles and zeros
            warning off all % avoid lots of ill-conditioning error messages
            nflin=100;      % number of frequency samples in linear region (high freq)
            alp=1.15;        % freq ratio increment in low freq region
            f0=10*2*pi/f;  % minimum interesting frequency (10 Hz in radians)
            fx=pi/nflin/(alp-1);    % boundary between log and linear portions
            if fx<=f0 || f0>=pi
                fif=linspace(0,pi,nflin);
            elseif fx>pi
                fif=[0 logspace(log10(f0),log10(pi),ceil(log10(pi/f0)/log10(alp)))];
            else
                nlin=ceil((pi-fx)*nflin/pi);
                fif=[0 logspace(log10(f0),log10(fx),ceil(log10(fx/f0)/log10(alp))) linspace(fx+(pi-fx)/nlin,pi,nlin-1)];
            end
            h0=abs(polyval(sb,1i*fif*f)./polyval(sa,1i*fif*f));   % target magnitude spectrum
            h0tol=max(h0)*1e-4;                                     % absolute gain tolerance
            % initialize with impulse invariance
            [bz,az]=impinvar(sb,sa,f);
            hj=freqz(bz,az,fif);
            maxj=max(abs(db(abs(hj)+h0tol)-db(h0+h0tol)));
            % or else with bilinear
            [ifb,ifa]=bilinear(sb,sa,f);
            hn=freqz(ifb,ifa,fif);
            maxi=max(abs(db(abs(hn+h0tol))-db(h0+h0tol)));
            if maxi<maxj || maxj>10 % accept bilinear if it is better or if imp inv is bad
                maxj=maxi;
                bz=ifb;
                az=ifa;
                hj=hn;
            end
            pat0=sb(end)==0;        % we have a zero at DC
            if pat0
                fif(1)=[];          % eliminate DC as a probe frequency
                hz1=1-exp(-1i*fif); % response of zero at z=1
                h0=h0(2:end)./abs(hz1); % remove a zero at z=1 from the target
                hj=hj(2:end)./hz1; % remove a zero at z=1 from the initial phase
            end
            upd=0;
            for mm=1:length(sa)     % maximum number of poles
                for nn=1:mm         % number of zeros is always less than number of poles
                    hn=hj;
                    j=0;
                    for i=1:30          % iterate up to 30 times (usually less)
                        h=h0.*exp(1i*angle(hn));
                        [ifb,ifa]=invfreqz(h,fif,nn,mm,[],10);
                        hn=freqz(ifb,ifa,fif);
                        maxi=max(abs(db(abs(hn+h0tol))-db(h0+h0tol)));
                        if maxi<maxj
                            maxj=maxi;
                            bz=ifb;
                            az=ifa;
                            hj=hn;
                            j=i;
                            upd=1;
                        end
                        if i>j+5    % quit if no improvement in last five iterations
                            break
                        end
                    end
                end
            end
            if upd
                bz=conv(bz,[1 -1]); % restore the zero at z=0
            end
            warning on all
            if si>0
                ient=ient+1; % cache entry number
                [jdum,jent]=min(fixz(:,1));     % find least recently used cache entry
                fixz(jent,1)=ient;  % flag as recently used
                fixz(jent,2)=f*nspty+si;   % save frequency/spectrum code
                baz{jent,1}=bz;
                baz{jent,2}=az;
            end
        end
    end
end
switch m1
    case 'z'
        b=bz;
        a=az;
    case 't'
        if nargin<5 || ~numel(zi)
            [b,a]=randfilt(bz,az,n);
        else
            [b,a]=randfilt(bz,az,n,zi);
        end
    case 'i'
        if nargin<5 || ~numel(zi)
            if nargin<4 || ~numel(n)  % determine n to include 1 - 1e-8 of the energy
                n=ceil(-fs*log(1e4)/max(real(sp)));
            end
            [b,a]=filter(bz,az,[1; zeros(n-1,1)]);
        else
            [b,a]=filter(bz,az,zeros(n,1),zi);
        end
    case 'm'
        b = abs(h);
    case 'f'
        b = h;
    case 'd'
        b = db(abs(h));
    case 'e'
        b=db(abs(h).*f*log(10)); % convert to power per decade in dB
    case 'l'
        b=h.*conj(h).*f*log(10); % convert to power per decade
    case 'p'
        b=h.*conj(h);
    case 's'
        b=sb;
        a=sa;
    otherwise
        error('Output format %s not implemented',m1);
end

% plot data
if ~nargout || ~strcmp(m,lower(m))
    if strcmp(m,lower(m)) % if no upper case letters
        m='MLT';
    end
    if ~any(m=='Q') && ~any(m=='A') && ~any(m=='E') && ~any(m=='F') && (m1~='t' || ~any(m=='W') && ~any(m=='S'))
        m(end+1)='M';  % default plot type
    end
    nfig=0;
    paz=any(m1=='itz');  % plot discrete time result
    pas=~paz || any(m=='T'); % plot continuous time result
    if any(m=='M') || any(m=='Q') || any(m=='E') || any(m=='E') % draw a frequency response plot
        clf;
        nfig=1;
        pam=any(m=='M');  % magnitude response
        pae=(any(m=='E') || any(m=='F')) && paz;  % magnitude response error
        paq=any(m=='Q');  % phase response
        pat=any(m=='T') && paz;  % include target spectrum
        pal=any(m=='L') || (~paz && ~any(m=='U'));  % log frequency axis
        if any(m1=='itz')
            fs=f;  % save the sample frequency
            apz=abs([sp sz]); % and determine the frequency range
            apz(apz==0)=[]; % ignore zero frequency poles/zeros
            if ~numel(apz)
                apz=[100 5000];
            elseif length(apz)==1
                apz=[apz/10 apz*10];
            end
            if pal
                f=logspace(log10(min([fs/1000 apz*0.05/pi])),log10(fs/2),200);
            else
                f=linspace(min([fs/1000 apz*0.05/pi]),fs/2,200);
            end
        end
        hs=freqs(sb,sa,2*pi*f);
        if paz
            hz=freqz(bz,az,f,fs);
        end

        axh=[];
        nax=pam+pae+paq; % number of axis sets
        titex='';
        if pam
            if nax>1
                subplot(nax,1,1);
            end
            if paz
                if pas
                    plot(f,db(abs(hs)),'--r',f,db(abs(hz)),'-b')
                    ymax=max(db(abs([hs hz])))+1;
                    ymin=min(db(abs([hs hz])))-1;
                    titex=' ( : = s, - = z)';
                    %                     legend('s-domain','z-domain','location','best');
                else
                    plot(f,db(abs(hz)),'-b')
                    ymax=max(db(abs(hz)))+1;
                    ymin=min(db(abs(hz)))-1;
                end
            else
                plot(f,db(abs(hs)),'-b')
                ymax=max(db(abs(hs)))+1;
                ymin=min(db(abs(hs)))-1;
            end
            if pal
                set(gca,'Xscale','log');
            end
            set(gca,'YLim',[max(ymin,ymax-60) ymax]);
            axisenlarge([-1 1.05]);
            xlabel(['Frequency (' xticksi 'Hz)']);
            ylabel('Gain (dB)');

            if si>0
                title(sprintf('Type %d: %s%s',si,spty{si},titex));
            end
            axh(end+1)=gca;
        end
        if pae
            if nax>1
                subplot(nax,1,1+pam);
            end
            if any(m=='F')
                dbflr=40;
                dbfl=max(abs(hs))*10^(-dbflr/20); % make a floor 40 dB below peak
                dbflt=sprintf(' (floor@-%.0fdB)',dbflr);
            else
                dbfl=0;
                dbflt='';
            end
            dberr=db(abs(hz)+dbfl)-db(abs(hs)+dbfl);
            plot([f(1) f(end)],[0 0],':k',f,dberr,'-b')
            if pal
                set(gca,'Xscale','log');
            end
            axisenlarge([-1 -1.05]);
            %             set(gca,'XLim',[min(f) max(f)]);
            xlabel(['Frequency (' xticksi 'Hz)']);
            ylabel('Gain Error (dB)');
            texthvc(0.95,0.95,sprintf('Err < %.1f dB%s',max(abs(dberr)),dbflt),'RTk');
            if si>0 && ~pam
                title(sprintf('Type %d: %s%s',si,spty{si},titex));
            end
            axh(end+1)=gca;
        end
        if paq
            if nax>1
                subplot(nax,1,nax);
            end
            if paz
                if pas
                    plot(f,angle(hs),'--r',f,angle(hz),'-b')
                else
                    plot(f,angle(hz),'-b')
                end
            else
                plot(f,angle(hs),'-b')
            end
            if pal
                set(gca,'Xscale','log');
            end
            axisenlarge([-1 -1.05]);
            %             set(gca,'XLim',[min(f) max(f)]);
            xlabel(['Frequency (' xticksi 'Hz)']);
            ylabel('Phase (rad)');
            if si>0 && nax==1
                title(spty{si});
            end
            axh(end+1)=gca;
        end
        if nax>1
            linkaxes(axh,'x');
        end
    end
    if any(m=='A') % plot complex plane
        if nfig
            figure();
        end
        clf;
        nfig=1;
        if pas
            if paz
                subplot(121);
            end
            plot(real(sp),imag(sp),'xb',real(sz),imag(sz),'ob');
            axis equal;
            xlim=get(gca,'xlim');
            xlim(1)=min(xlim(1),-1000);
            xlim(2)=max(xlim(2),1000);
            ylim=get(gca,'ylim');
            ylim(1)=min(ylim(1),-1000);
            ylim(2)=max(ylim(2),1000);
            axis([xlim ylim]);
            hold on
            plot(xlim,[0 0],':r',[0 0],ylim,':r');
            hold off
            title('s-plane');
        end
        if paz
            if pas
                subplot(122);
            end
            axl=max(abs([1.1; az(:);bz(:)]));
            t=linspace(0,2*pi);
            rtzb=roots(bz);
            rtza=roots(az);
            plot(cos(t),sin(t),':r',[-1 0; 1 0],[0 -1; 0 1],':r',real(rtza),imag(rtza),'xb',real(rtzb),imag(rtzb),'ob');
            axis equal;
            axis([-1 1 -1 1]*axl);
            title('z-plane');
        end
    end
    if any(m=='W') && any(m1=='it') % plot waveform
        if nfig
            figure();
        end
        clf;
        nfig=1;
        plot((1:length(b))/fs,b,'-b');
        xlabel(['Time (' xticksi 's)']);
        if si>0
            title(spty{si});
        end
    end
    if any(m=='S') && any(m1=='it') && numel(b)>0.1*fs % plot spectrogram
        if nfig
            figure();
        end
        clf;
        if any(m=='L')
            sm='pJcwl';
        else
            sm='pJcw';
        end
        spgrambw(b,fs,sm);
        if si>0
            title(spty{si});
        end
    end
end