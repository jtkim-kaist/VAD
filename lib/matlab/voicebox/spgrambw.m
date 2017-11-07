function [t,f,b]=spgrambw(s,fs,varargin)
%SPGRAMBW Draw spectrogram [T,F,B]=(s,fs,mode,bw,fmax,db,tinc,ann)
%
%  Usage: (1) spgrambw(s,fs,'pJcw')     % Plot spectrogram with my favourite set of options
%
%         (2) spgrambw(s,fs,'PJcwm',50,[100 2000])    % Plot narrow-band spectrogram on mel scale
%                                                       from 100 to 2000 mel in power/mel units
%
%         (3)     [t,f,b]=spgrambw(s,fs,'p');        % calculate the spectrogram without plotting
%                 imagesc(t,f,10*log10(b'));         % plot it manually
%                 axis('xy');
%
%         (4) ninc=0.0045*fs;           % Frame increment for BW=200 Hz (in samples)
%             nwin=2*ninc;              % Frame length (in samples)
%             win=hamming(nwin);        % Analysis window
%             k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
%             sf=abs(rfft(enframe(s,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array                
%             spgrambw(sf,[fs/ninc 0.5*(nwin+1)/fs fs/nwin],'Jc',bw);  % Plot spectrum array
%
%         For examples of the many options available see:
%         http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/tutorial/spgrambw/spgram_tut.pdf
%
%  Inputs:  S         speech signal, or single-sided power spectrum array, S(NT,NF), in power per Hz
%           FS        sample fequency (Hz) or [FS T1] where T1 is the time of the first sample
%                     or, if s is a matrix, [FS T1 FINC F1] where FS is the frame rate, T1 is
%                     the time of the first sample, FINC is the frequency increment and F1 the
%                     frequency of the first column.
%           MODE      optional character string specifying options (see list below)
%           BW        bandwidth resolution in Hz (DFT window length = 1.81/BW)[default: 200]
%           FMAX      frequency range [Fmin Fstep Fmax]. If Fstep is omitted
%                     it is taken to be (Fmax-Fmin)/257, if Fmin is also omitted it is taken
%                     to be 0 (or 20Hz for mode l), if all three are omitted Fmax is taken to be FS/2.
%                     If modes m, b, e or l are specified then the units are in mel, bark or erb or
%                     log10(Hz); this can be over-ridden by the 'h' option.
%           DB        either dB-range relative to peak or [dB-min dB-max] for plotting (and B output if 'D' given [default: 40]
%           TINC      output frame increment in seconds [0 or missing uses default=0.45/BW]
%                     or [TFIRST TLAST] or [TFIRST TINC TLAST] where TFIRST/TLAST are the times
%                     of first/last frames
%           ANN       annotation cell array: each row contains either
%                     {time 'text-string' 'font'} or {[t_start t_end] 'text-string' 'font'} where
%                     the time value is in seconds with s(n) at time offset+n/fs. The font column can
%                     omitted in which case the system font will be used. MATLAB cannot cope with
%                     unicode so I recommend the SILDoulosIPA (serifed) or SILSophiaIPA (sans) fonts
%                     for phonetic symbols; these are now a little hard to find.
%
% Outputs:  T(NT)        time axis values (in seconds). Input sample s(n) is at time offset+n/fs.
%           F(NF)        frequency axis values in Hz or, unless mode=H, other selected frequency units
%                        according to mode: m=mel, l=log10(Hz), b=bark,e=erb-rate
%           B(NT,NF)     spectrogram values in power per x where x depends on the 'pPmbel' options
%                        clipped to DB range if 'D' option and in dB if 'd' option.
%
% MODE:  'p' = output power per decade rather than power per Hz [preemphasis]
%        'P' = output power per mel/bark/erb according to y axis scaling
%        'd' = output B array is in dB rather than power
%        'D' = clip the output B array to the limits specified by the "db" input
%
%        'n' = use nearest frequency bin instead of interpolating
%
%        'm' = mel scale
%        'b' = bark scale
%        'e' = erb scale
%        'l' = log10 Hz frequency scale
%        'f' = label frequency axis in Hz rather than mel/bark/...
%
%        'h' = units of the FMAX input are in Hz instead of mel/bark
%              [in this case, the Fstep parameter is used only to determine
%               the number of filters]
%        'H' = express the F output in Hz instead of mel/bark/...
%
%        'g' = draw a graph even if output arguments are present
%        'j' = jet colourmap
%        'J' = "thermal" colourmap that is linear in grayscale. Based on Oliver Woodford's
%                 real2rgb at http://www.mathworks.com/matlabcentral/fileexchange/23342
%        'i' = inverted colourmap (white background)
%        'c' = include a colourbar as an intensity scale
%        'w' = draw the speech waveform above the spectrogram
%        'a' = centre-align annotations rather than left-aligning them
%        't' = add time markers with annotations
%
% The BW input gives the 6dB bandwidth of the Hamming window used in the analysis.
% Equal amplitude frequency components are guaranteed to give separate peaks if they
% are this far apart. This value also determines the time resolution: the window length is
% 1.81/BW and the low-pass filter applied to amplitude modulations has a 6-dB bandwidth of
% BW/2 Hz.
%
% The units are power per Hz unless the 'P' option is given in which case power
% per displayed unit is used or power per decade for the 'l' option.

%%%% BUGS %%%%%%
% * allow ANN rows to be a mixture of intervals and instants
% * allow multiple ANN rows
% * frequency axis labels ofter start at -0 instead of 0.
% * Do not use triangular interpolation if the output frequencies are the same as an FFT
% * Place as many subticks as will fit beyond the last tick with the 'f' option
% * Use a special subtick pattern between ticks that are powers of 10 using the 'f' option
% * Future options:
%       ['q' = constant q transform]
%       ['k' = add a piano keyboard to the frequency scale]
%       ['z' = use a bipolar colourmap for a matrix input with negative values]

%      Copyright (C) Mike Brookes 1997-2011
%      Version: $Id: spgrambw.m 6803 2015-09-12 09:31:44Z dmb $
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
persistent tcmap
if isempty(tcmap)
    % modified thermal with better grayscale linearity
    tcmap=[ 0 0 0; 7 0 17; 14 0 33; 21 0 50; 29 0 67; 36 0 84; 43 0 100; 50 0 117;
        57 0 134; 64 0 150; 72 0 167; 80 3 164; 89 7 156; 97 11 149; 106 15 142; 114 19 134;
        123 23 127; 131 27 119; 140 31 112; 149 35 105; 157 39 97; 166 43 90; 174 47 82;
        183 51 75; 192 55 68; 200 59 60; 209 63 53; 217 67 45; 226 71 38; 234 75 31;
        243 79 23; 252 83 16; 255 88 12; 255 95 12; 255 102 11; 255 109 11; 255 116 10;
        255 123 10; 255 130 9; 255 137 9; 255 144 8; 255 151 8; 255 158 7; 255 165 7;
        255 172 6; 255 179 6; 255 186 5; 255 193 4; 255 200 4; 255 207 3; 255 214 3; 255 221 2;
        255 228 2; 255 235 1; 255 242 1; 255 249 0; 255 252 22; 255 252 55; 255 253 88;
        255 253 122; 255 254 155; 255 254 188; 255 255 222; 255 255 255]/255;
end
if nargin<2
    error('Usage: SPGRAMBW(s,fs,mode,bw,fmax,db,tinc)');
end
%SPGRAMBW Draw grey-scale spectrogram [T,F,B]=(s,fs,mode,bw,fmax,db,tinc)
%
% first decode the input arguments
%
if size(s,1)==1
    s=s(:);   % force to be a column vector (unless it is a matrix)
end
[ns1,ns2]=size(s);
ap=zeros(1,6);
j=2;
if numel(fs)<2
    fs(2)=1/fs(1);  % first sample or frame is at time 1/fs
end
for i=1:length(varargin)
    if ischar(varargin{i})
        ap(1)=i;
    else
        ap(j)=i;
        j=j+1;
    end
end
if ap(1) && ~isempty(varargin{ap(1)})
    mode=varargin{ap(1)};
else
    mode='';  % default mode
end
if ap(2) && ~isempty(varargin{ap(2)})
    bw=varargin{ap(2)};
else
    bw=200;
end
if ap(3) && ~isempty(varargin{ap(3)})
    fmax=varargin{ap(3)};
else
    fmax=[];
end
if ap(4) && ~isempty(varargin{ap(4)})
    db=varargin{ap(4)};
else
    db=40;
end
if ap(5) && ~isempty(varargin{ap(5)})
    tinc=varargin{ap(5)};
else
    tinc=0;
end
switch numel(tinc)
    case 1
        tinc=[tinc -Inf Inf];
    case 2
        tinc=[0 tinc];
    otherwise
        tinc=tinc([2 1 3]);
end
if tinc(1)<=0
    tinc(1)=1.81/(4*bw); % default frame increment
end
if ap(6)
    ann=varargin{ap(6)};
else
    ann=[];
end

% now sort out the mode flags

mdsw='  ';           % [yscale preemph]
for i=1:length(mode)
    switch mode(i)
        case {'l','m','b','e'}
            mdsw(1)=mode(i);
        case {'p','P'}
            mdsw(2)=mode(i);
    end
end
if mdsw(2)=='P'
    mdsw(2)=mdsw(1);        % preemphasis is scaling dependent
end
%
% sort out the frequency axis
%
flmin=30;                   % min frequency for 'l' option
nfrq=257;                   % default number of frequency bins
if ns2==1
    fnyq=fs(1)/2;           % default upper frequency limit is fs/2
else                        % input is a power spectrum
    if numel(fs)<3
        fs(3)=fs(1)*0.25;   % default increment is 0.25 times frame increment
    end
    if numel(fs)<4
        fs(4)=0;            % first freq bin is DC by default
    end
    fnyq=fs(4)+(ns2-1)*fs(3);  % default upper frequency limit is highest supplied frequency
end

if ~numel(fmax)             % no explicit frequency range
    switch mdsw(1)
        case 'l'
            fx=linspace(log10(flmin),log10(fnyq),nfrq);   % 20  Hz to Nyquist
        case 'm'
            fx=linspace(0,frq2mel(fnyq),nfrq);   % DC to Nyquist
        case 'b'
            fx=linspace(0,frq2bark(fnyq),nfrq);   % DC to Nyquist
        case 'e'
            fx=linspace(0,frq2erb(fnyq),nfrq);   % DC to Nyquist
        otherwise   % linear Hz scale
            fx=(0:nfrq-1)*fnyq/(nfrq-1);
    end
else
    if any(mode=='h')
        switch mdsw(1)
            case 'l'
                fmaxu=log10(fmax);   % 20  Hz to Nyquist
            case 'm'
                fmaxu=frq2mel(fmax);   % DC to Nyquist
            case 'b'
                fmaxu=frq2bark(fmax);   % DC to Nyquist
            case 'e'
                fmaxu=frq2erb(fmax);   % DC to Nyquist
            otherwise
                fmaxu=fmax;  % linear Hz scale
        end
    else
        fmaxu=fmax;                 % already in the correct units
    end
    if numel(fmax)<2   % only max value specified
        if mdsw(1)=='l'
            fx=linspace(log10(flmin),fmaxu,nfrq);   % 20  Hz to fmax
        else
            fx=linspace(0,fmaxu,nfrq);   % DC to fmax
        end
    elseif numel(fmax)<3 % min and max values specified
        fx=linspace(fmaxu(1),fmaxu(2),nfrq);   % fmin to fmax
    else
        fmaxu(2)=fmax(2)*(fmaxu(3)-fmaxu(1))/(fmax(3)-fmax(1)); % scale the step size appropriately
        fx=fmaxu(1):fmaxu(2):fmaxu(3);   % fmin to fmax in steps of finc
        nfrq=length(fx);
    end
end
switch mdsw(1)          % convert the frequency range to Hz
    case 'l'
        f=10.^fx;
        frlab='log_{10}Hz';
        frlabf='log';
        frq2y=@log10;
        y2frq=@(x) 10.^x;
    case 'm'
        f=mel2frq(fx);
        frlab='Mel';
        frlabf='Mel';
        frq2y=@frq2mel;
        y2frq=@mel2frq;
    case 'b'
        f=bark2frq(fx);
        frlab='Bark';
        frlabf='Bark';
        frq2y=@frq2bark;
        y2frq=@bark2frq;
    case 'e'
        f=erb2frq(fx);
        frlab='Erb-rate';
        frlabf='Erb';
        frq2y=@frq2erb;
        y2frq=@erb2frq;
    otherwise
        f=fx;
        frlab='Hz';
        frq2y=@(x) x;
        y2frq=@(x) x;
end
if ~any(mode=='H')
    f=fx;               % give output frequencies in native units instead of Hz unless 'H' is specified
end
%
% now calculate the spectrogram
%
if ns2==1   % input is a speech signal vector
    winlen = fix(1.81*fs(1)/bw);   % window length
    win=0.54+0.46*cos((1-winlen:2:winlen)*pi/winlen);  % Hamming window
    ninc=max(round(tinc(1)*fs(1)),1);                 % window increment in samples
    %  we need to take account of minimum freq increment + make it exact if possible
    fftlen=pow2(nextpow2(4*winlen));        % enough oversampling to get good interpolation
    win=win/sqrt(sum(win.^2));              % ensure window squared sums to unity
    ix1=max(round((tinc(2)-fs(2))*fs(1)-(winlen-3)/2),1); % first sample required
    ix2=min(ceil((tinc(3)-fs(2))*fs(1)+(winlen+1)/2),ns1); % last sample required
    [sf,t]=enframe(s(ix1:ix2),win,ninc);
    t=fs(2)+(t+ix1-2)/fs(1);                         % time axis
    b=rfft(sf,fftlen,2);
    b=b.*conj(b)*2/fs(1);          % Power per Hz
    b(:,1)=b(:,1)*0.5;   % correct for no negative zero frequency to double the power
    b(:,end)=b(:,end)*0.5;   % correct for no negative nyquist frequency to double the power
    fb=(0:fftlen/2)*fs(1)/fftlen; % fft bin frequencies
    fftfs=fs(1);
else
    b=s;
    t=fs(2)+(0:ns1-1)/fs(1);  % frame times
    fb=fs(4)+(0:ns2-1)*fs(3);
    fftlen=[ns2 fs(3) fs(4)]; % for filtbankm: ns2=# input freq bins, freq increment (Hz), first bin freq (Hz)
    fftfs=0;
    %     fftlen=2*(ns2-1);  % assume an even length fft
    %     fftfs=fftlen*fs(3);
end
nfr=numel(t);                   % number of frames
dblab='Power/Hz';
switch mdsw(2)
    case {'p','l'}
        b=b.*repmat(fb*log(10),nfr,1);       % convert to power per decade
        dblab='Power/Decade';
    case 'm'
        b=b.*repmat((700+fb)*log(1+1000/700)/1000,nfr,1);       % convert to power per mel
        dblab='Power/Mel';
    case 'b'
        b=b.*repmat((1960+fb).^2/52547.6,nfr,1);       % convert to power per bark
        dblab='Power/Bark';
    case 'e'
        b=b.*repmat(6.23*fb.^2 + 93.39*fb + 28.52,nfr,1);       % convert to power per erb
        dblab='Power/Erb-rate';
end
%
% Now map onto the desired frequency scale
%
if any(mode=='n')
    fbopt=['cushn' mdsw(1)];
else
    fbopt=['cush' mdsw(1)];
end
b=b*filtbankm(nfrq,fftlen,fftfs,fx(1),fx(end),fbopt)';

if ~nargout || any(mode=='g') ||  any(lower(mode)=='d')
    if numel(db)<2          % find clipping limits
        plim=max(b(:))*[0.1^(0.1*db) 1];
    else
        plim=10.^(0.1*db(1:2));
    end
    if plim(2)<=0
        plim(2)=1;
    end
    if plim(1)<=0 || plim(1)==plim(2)
        plim(1)=0.1*plim(2);
    end
    if ~nargout || any(mode=='g')
        bd=10*log10(max(b,max(b(:)*1e-30)));  % save an unclipped log version for plotting
    end
    if any(mode=='D')
        b=min(max(b,plim(1)),plim(2)); % clip the output
    end
    if any(mode=='d')
        b=10*log10(b);    % output the dB version
    end
end
% now plot things
if ~nargout || any(mode=='g')
    cla;  % clear current axis
    imagesc(t,fx,bd');
    axis('xy');
    set(gca,'tickdir','out','clim',10*log10(plim));
    if any(mode=='j')
        colormap('jet');
        map=colormap;
    elseif any(mode=='J')
        map=tcmap;
    else
        map = repmat((0:63)'/63,1,3);
    end
    if any(mode=='i')               % 'i' option = invert the colourmap
        map=map(64:-1:1,:);
    end
    colormap(map);
    if any(mode=='c')                % 'c' option = show a colourbar
        colorbar;
        cblabel([dblab ' (dB)']);
    end
    %
    % Now check if annotations or a waveform are required
    %
    dotaw=[((any(mode=='t') && size(ann,2)>1) || size(ann,2)==1) size(ann,2)>1 (any(mode=='w') && ns2==1)];
    ylim=get(gca,'ylim');
    if  any(dotaw)
        yrange = ylim(2)-ylim(1);
        zlim=ylim;
        toptaw=cumsum([0 dotaw.*[0.05 0.05 0.1]]*yrange)+ylim(2);
        zlim(2)=toptaw(4);
        set(gca,'ylim',zlim,'color',map(1,:));
        if dotaw(3)        % Plot the waveform
            six=min(max(floor((get(gca,'xlim')-fs(2))*fs(1))+[1 2],1),ns1);
            smax=max(s(six(1):six(2)));
            smin=min(s(six(1):six(2)));
            if smax==smin
                smax=smax+1;
                smin=smin-1;
            end
            srange=smax-smin;
            hold on
            plot(fs(2)+(six(1)-1:six(2)-1)/fs(1),(s(six(1):six(2))-smin)/srange*0.9*(toptaw(4)-toptaw(3))+toptaw(3),'color',map(48,:))
            hold off
        end
        if dotaw(1) || dotaw(2)
            tmk=cell2mat(ann(:,1));
            tmksel=tmk(:,1)<=t(end) & tmk(:,end)>=t(1);
            yix=1+[tmk(tmksel,1)<t(1) ones(sum(tmksel),2) tmk(tmksel,end)>t(end)]';
            tmk(:,1)=max(tmk(:,1),t(1));  % clip to axis limits
            tmk(:,end)=min(tmk(:,end),t(end));
        end
        if dotaw(1) && any(tmksel)  % draw time markers
            ymk=toptaw(1:2)*[0.8 0.4;0.2 0.6];
            switch size(tmk,2)
                case 0
                case 1      % isolated marks
                    hold on
                    plot([tmk(tmksel) tmk(tmksel)]',repmat(ymk',1,sum(tmksel)),'color',map(48,:));
                    hold off
                otherwise % draw durations

                    hold on
                    plot(tmk(tmksel,[1 1 2 2])',ymk(yix),'color',map(48,:));
                    hold off
            end
        end
        if dotaw(2) && any(tmksel) % print annotations
            if any(mode=='a')
                horal='center';
                tmk=(tmk(:,1)+tmk(:,end))*0.5;
            else
                horal='left';
                tmk=tmk(:,1);
            end
            if size(ann,2)>2
                font='Arial';
                for i=1:size(ann,1)
                    if tmksel(i)
                        if ~isempty(ann{i,3})
                            font = ann{i,3};
                        end
                        text(tmk(i),toptaw(2),ann{i,2},'color',map(48,:),'fontname',font,'VerticalAlignment','baseline','HorizontalAlignment',horal);
                    end
                end
            else
                for i=1:size(ann,1)
                    if tmksel(i)
                        text(tmk(i),toptaw(2),ann{i,2},'color',map(48,:),'VerticalAlignment','baseline','HorizontalAlignment',horal);
                    end
                end
            end
        end
    end
    xlabel(['Time (' xticksi 's)']);
    if any(mode=='f') && ~strcmp(frlab,'Hz')
        ylabel([frlabf '-scaled frequency (Hz)']);
        ytickhz(frq2y,y2frq);
    else
        ylabel(['Frequency (' yticksi frlab ')']);
    end
    ytick=get(gca,'YTick');
    ytickl=get(gca,'YTickLabel');
    msk=ytick<=ylim(2);
    if any(~msk)
        set(gca,'YTick',ytick(msk),'YTickLabel',ytickl(msk));
    end
end

function ytickhz(frq2y,y2frq)
% label non linear y frequency axis
%
% Bugs/Suggestions:
% * Add a penalty for large numbers (e.g. 94 is less "round" than 11)
% * possibly add subticks at 1:2:5 if boundaries are 1 and 10
% * could treat subtick allocation specially if bounding lables are both powers of 10
%   and work in log spacing rather than spacing directly

% algorithm constants

seps=[0.4 1 3 6]; % spacings: (a) min subtick, (b) min tick, (c) min good tick, (d) max good tick
ww=[0.5 0.6 0.8 0.1 0.3 0.3 0.2];  % weight for (a) last digit=5, (b) power of 10, (c) power of 1000, (d) equal spacing, (e) 1:2:5 labels (f) <seps(3) (g) >seps(4)
nbest=10; % number of possibilities to track

prefix={'y','z','a','f','p','n','µ','m','','k','M','G','T','P','E','Z','Y'};

ah=gca;
getgca=get(ah);  % Get original axis properties
set(ah,'Units','points','FontUnits','points');
getgcac=get(ah);  % Get axis properties in points units
set(ah,'Units',getgca.Units,'FontUnits',getgca.FontUnits); % return to original values
ylim=getgca.YLim;
yrange=ylim*[-1;1];
chsz= yrange*getgcac.FontSize/getgcac.Position(4); % char height in Y-units
% divide the y-axis up into bins containing at most one label each
maxl=ceil(2*yrange/chsz);  % max number of labels

% candidate array [cand(:,[1 2])/1000 cand(:,5) cand(:,6)/1000 cand(:,[7 8])]
% 1,2=y limits, 3,4=log limits, 5=Hz, 6=cost, 7=mantissa, 8=exponent, 9=sig digits, 10=y-position
cand=zeros(maxl+2,10);
yinc=(yrange+chsz*0.0002)/maxl;  % bin spacing (allowing for a tiny bit to ensure the ends are included)
cand(2:end-1,2)=ylim(1)+yinc*(1:maxl)'-chsz*0.0001;
cand(3:end-1,1)=cand(2:end-2,2);
cand(2,1)=cand(2,2)-yinc;
cand(2:end-1,1:2)=y2frq(max(cand(2:end-1,1:2),0));

% find the "roundest" number in each interval
% first deal with intervals containing zero
cand([1 maxl+2],6)=-1;
cand(2,9)=(cand(2,1)<=0);  % mask out interval contaiing zero
cand(2,6)=-cand(2,9);
msk=cand(:,6)==0;  % find rows without a cost yet
cand(msk,3:4)=log10(cand(msk,1:2));
% find powers of 1000
loglim=ceil(cand(:,3:4)/3);
msk=loglim(:,2)>loglim(:,1);
if any(msk)
    xp=loglim(msk,1);
    wuns=ones(length(xp),1);
    cand(msk,5:9)=[1000.^xp wuns-ww(3) wuns 3*xp wuns];
end
% find powers of 10
loglim=ceil(cand(:,3:4));
msk=~msk & (loglim(:,2)>loglim(:,1));
if any(msk)
    xp=loglim(msk,1);
    wuns=ones(length(xp),1);
    cand(msk,5:9)=[10.^xp wuns-ww(2) wuns xp wuns];
end
% find value with fewest digits
msk=cand(:,6)==0;  % find rows without a cost yet
maxsig=1-floor(log10(10^min(cand(msk,3:4)*[-1;1])-1)); % maximum number of significant figures to consider
pten=10.^(0:maxsig-1);   % row vector of powers of ten
noten=10.^(-floor(cand(msk,3))); % exponent of floating point representation of lower bound
sigdig=sum((ceil(cand(msk,2).*noten*pten)-ceil(cand(msk,1).*noten*pten))==0,2); % number of digits common to the interval bounds
lowman=ceil(cand(msk,1).*noten.*10.^sigdig);
midman=10*floor(lowman/10)+5;
highman=ceil(cand(msk,2).*noten.*10.^sigdig);
mskman=midman>=lowman & midman<highman;   % check if we can include a manitssa ending in 5
lowman(mskman)=midman(mskman);
cand(msk,6:9)=[sigdig+1 lowman floor(cand(msk,3))-sigdig sigdig+1];
cand(msk,5)=cand(msk,7).*10.^cand(msk,8);
cand(msk,6)=cand(msk,6)-(mod(cand(msk,7),10)==5)*ww(1);
cand(2:end-1,10)=frq2y(cand(2:end-1,5));
cand([1 maxl+2],10)=ylim + seps(4)*chsz*[-1 1]; % put imaginary labels at the optimum spacing beyond the axes
% [cand(:,[1 2 5])/1000 cand(:,[6 7 8 9])]

% Now do n-best DP to find the best sequence

ratint=[8/5 25/10 0 0 4/3];
costs=Inf(nbest,maxl+2); % cumulative path costs
costs(1,1)=0; % starting node only has one option
prev=ones(nbest,maxl+2); % previous label in path
labcnt=zeros(nbest,maxl+2); % number of labels in path
for i=2:maxl+2
    ntry=nbest*(i-1); % number of previous options
    prevc=reshape(repmat(1:i-1,nbest,1),ntry,1); % previous candidate
    prevprev=1+floor((prev(1:ntry)'-1)/nbest); % previous previous candidate
    msk=prevprev>1+(maxl+2)*(i==maxl+2); % mask for label triplets
    labcnti=labcnt(1:ntry)+1;
    disti=(cand(i,10)-cand(prevc,10))/chsz; % distance to previous label in characters
    costa=max(seps(3)-disti,0)*ww(6)+max(disti-seps(4),0)*ww(7);
    incri=(cand(i,5)-cand(prevc,5)); % label increment
    incrj=(cand(i,5)-cand(prevprev,5)); % double label increment
    if any(msk)
        costa(msk)=costa(msk)- ww(4)*(abs(incrj(msk)-2*incri(msk))<0.01*incri(msk));
        if cand(i,7)==1 || cand(i,7)==2 || cand(i,7)==5 % look for labels 1:2:5
            costa(msk)=costa(msk)- ww(5)*(abs(incrj(msk)-ratint(cand(i,7))*incri(msk))<0.01*incri(msk));
        end
    end
    costa(disti<seps(2))=Inf;
    costi=(costs(1:ntry).*max(labcnt(1:ntry),1)+costa'+cand(i,6))./labcnti;
    [sc,isc]=sort(costi);
    isc=isc(1:nbest);
    costs(:,i)=sc(1:nbest)';
    prev(:,i)=isc';
    labcnt(:,i)=labcnti(isc)';
end

% now traceback the best sequence

% fprintf('Traceback\n\n');
ichoose=0;
labchoose=[];
for i=1:nbest
    if labcnt(i,maxl+2)>1 && costs(i,maxl+2)<Inf
        lablist=zeros(labcnt(i,maxl+2)-1,1);
        k=prev(i,maxl+2);
        for j=labcnt(i,maxl+2)-1:-1:1
            lablist(j)=1+floor((k-1)/nbest);
            k=prev(k);
        end
        %         fprintf('Cost=%8.2f :',costs(i,maxl+2));
        %         fprintf(' %g',cand(lablist,5))
        %         fprintf('\n');
        if ~ichoose || labcnt(ichoose,maxl+2)==1
            ichoose=i;
            labchoose=lablist;
        end
    end
end

% now create the labels

ntick=length(labchoose);
% sort out the subticks
subpos=[];
if ntick>=2
    for i=1:ntick-1
        clj=cand(labchoose(i:i+1),:);
        sprec=min(clj(1,8)+100*(clj(1,7)==0),clj(2,8)); % subtick precision
        spos=(clj(1,7)*10^(clj(1,8)-sprec):clj(2,7)*10^(clj(2,8)-sprec))*10^sprec;
        nsub=length(spos);
        if nsub==2
            spos=spos*[1 0.5 0;0 0.5 1];
            nsub=3;
        end
        if nsub>=3
            yspos=frq2y(spos);
            for kk=1:3 % try various subdivisions: every 1, 2 or 5
                k=kk+2*(kk==3);  % 1, 2 and 5
                if 2*k<=nsub-1 && ~mod(nsub-1,k)  % must divide exactly into nsub
                    if all((yspos(1+k:k:nsub)-yspos(1:k:nsub-k))>=(seps(1)*chsz)) % check they all fit in
                        subpos=[subpos yspos(1+k:k:nsub-k)];
                        if i==1
                            spos=(ceil(cand(2,1)/10^sprec):clj(1,7)*10^(clj(1,8)-sprec))*10^sprec;
                            nsub=length(spos);
                            yspos=frq2y(spos);
                            if nsub>=k+1 && all((yspos(nsub:-k:1+k)-yspos(nsub-k:-k:1))>=(seps(1)*chsz))
                                subpos=[subpos yspos(nsub-k:-k:1)];
                            end
                        elseif i==ntick-1
                            spos=(clj(2,7)*10^(clj(2,8)-sprec):floor(cand(end-1,2)/10^sprec))*10^sprec;
                            nsub=length(spos);
                            yspos=frq2y(spos);
                            if nsub>=k+1 && all((yspos(1+k:k:nsub)-yspos(1:k:nsub-k))>=(seps(1)*chsz))
                                subpos=[subpos yspos(1+k:k:nsub)];
                            end
                        end
                        break;
                    end
                end
            end
        end
    end
end
nsub=length(subpos);
tickpos=[cand(labchoose,10); subpos'];
ticklab=cell(ntick+nsub,1);
sipref=min(max(floor((sum(cand(labchoose,8:9),2)-1)/3),-8),8);
nzadd=cand(labchoose,8)-3*sipref;  % trailing zeros to add
digzer=cand(labchoose,7).*10.^max(nzadd,0); % label digits including trailing zeros
ndleft=cand(labchoose,9)+nzadd; % digits to the left of the decimal point
for i=1:ntick
    tickint=num2str(digzer(i));
    if nzadd(i)<0
        tickint=[tickint(1:ndleft(i)) '.' tickint(1+ndleft(i):end)];
    end
    ticklab{i} = sprintf('%s%s',tickint,prefix{sipref(i)+9});
end
for i=ntick+1:ntick+nsub
    ticklab{i}='';
end
[tickpos,ix]=sort(tickpos);
ticklab=ticklab(ix);

set(ah,'YTick',tickpos','YTickLabel',ticklab);

