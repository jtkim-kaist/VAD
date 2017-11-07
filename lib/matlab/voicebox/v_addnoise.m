function [z,p,fso]=v_addnoise(s,fsx,snr,m,nb,fsa)
%V_ADDNOISE Add noise at a chosen SNR [z,p,fso]=(s,fsx,snr,m,nb,fsa)
%
% Usage: (1) z=v_addnoise(s,fs,snr); % add white noise using P.56 for speech power with total ouput level at 0 dB
%        (2) z=v_addnoise(s,fs,snr,'',n,fn); % add noise from column vector n() at sample frequency fn with random start
%                                            % sample and wrapping around as required with a vorbis cross-fade
%        (3) z=v_addnoise(s,fs,snr,'AD');    % use A-weighting when calculating the SNR which is specified as a power ratio
%        (4) z=v_addnoise(s,fs,snr,'f',b,a); % generate noise using filter(b,a,randn(*,1)) but avoiding startup transient
%        (5) z=v_addnoise(s,fs,snr,'g',5);    % add speech-shaped noise (noise 5 from stdspectrum.m)and plot a graph
%
% Inputs:  s    input signal (column vector)
%        fsx    sampling frequency (Hz) or fso output from previous call
%        snr    target SNR in dB (or power ratio if 'D' option given) [default: Inf]
%          m    mode string [default = 'dxopEk']
%                (1) 'A'  use A-weighting when calculating SNR
%                (2) 'd'  SNR input and p output given in dB [default]
%                    'D'  SNR input and p output given as power ratio
%                (3) 'S'  no noise wrapping
%                    'w'  wrap noise if it is too short
%                    'W'  wrap noise with vorbis cross-fade of +-50 ms [default]
%                (4) 'b'  Go fom the beginning of the noise signal
%                    'o'  Add a random noise offset even if it means extra wraps [default]
%                    'O'  Add a random noise offset but minimize number of wraps
%                (5) 'p'  Use P.56 to calculate signal level [default]
%                    'e'  Use energy to calculate signal level
%                    'P'  Use P.56 to calculate noise level
%                    'E'  Use energy to calculate noise level [default]
%                (6) 'z'  Assume signal is already at 0 dB
%                    'Z'  Assume noise is already at 0 dB
%                (7) 'n'  make signal level = 0dB
%                    'N'  make noise level = 0dB
%                    't'  make sum of speech and noise levels = 0dB [default]
%                    'k'  preserve original signal power
%                    'K'  preserve original noise power
%                (8) 'x'  the output z contains input and noise as separate columns
%                (8) 'f'  Inputs nb and fsa specify a z-domain filter nb(z)/fsa(z) apllied to Gaussian noise
%                (9) 'g'  plot graph [default if no output arguments]
%         nb    noise signal or stdspectrum type or numerator of noise filter if 'f' option (default: white noise)
%        fsa    noise sample frequency [default fsx] or denominator of noise filter if 'f' option
%
% Outputs: z    noisy signal (single column unless 'x' option given
%          p    levels in dB or power (see 'dD' options): [s-in n-in s-out n-out]
%        fso    state output for fsx input to future call to allow s to be processed in blocks
%
% The noise can have a different sample rate from the signal (fsa and fsx
% inputs respectively) but if you will be re-using the noise multiple times
% it is more efficient to use v_resample() on the noise signal. This is the
% same as resample() but discards the initial and final filter transients.
% By using the fso output as the fsx input in a subsequent call, you can
% process a long speech signal in chuncks with the same noise. Note that
% the speech and noise levels determined fom the first chunck are used for
% all subsequent chuncks.

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: v_addnoise.m 5119 2014-09-11 07:22:12Z dmb $
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

ns=numel(s);        % number of speech samples
if isstruct(fsx)
    fs=fsx.fs;      % sample frequency
    nb=fsx.nb;        % noise samples or filter numerator
    genf=fsx.genf;  % filtered white noise flag
    nsc=fsx.nsc;    % cumulative sample count
    gm=fsx.gm; % gains for speech and noise
    p=fsx.p; % output powers
    if genf
        zof=fsx.zof;  % noise generating filter state
        fsa=fsx.fsa;  % filter denominator
    else
        ioff=fsx.ioff;  % offset in noise signal
    end
    %% decode input arguments
else
    fs=fsx;         % sample frequency
    nsc=0;          % cumulative sample count
    if nargin<3 || ~numel(snr)
        snr=Inf;
    end
    if nargin<4 || ~numel(m)
        m='';
    end
    if nargin<5 || ~numel(nb)
        nb=1;   % default is white noise
    end
    if any(m=='A') && (~any(m=='z') || ~any(m=='Z'))
        [awb,awa]=stdspectrum(2,'z',fs); % create an A-weighting filter
    end
    if any(m=='z')
        se=1;                       % speech power given as 0 dB
    elseif any(m=='A')
        sf=filter(awb,awa,s);        % apply A-weighting
        if any(m=='e')
            se=mean(sf.^2);         % use normal energy
        else
            se=activlev(sf,fs);     % or else P.56
        end
    else
        if any(m=='e')
            se=mean(s.^2);          % use normal energy
        else
            se=activlev(s,fs);      % or else P.56
        end
    end
    %% Now sort out the noise
    genf=any(m=='f') || ischar(nb) || numel(nb)==1;     % generate noise locally
    if genf   % generate noise locally
        if ~any(m=='f') % specify standard spectrum
            [nb,fsa]=stdspectrum(nb,'z',fs);
        end
        [dum1,zof,dum2,ne]=randfilt(nb,fsa,0);
        if any(m=='A')                                  % convolve with A-weighting to find power
            [dum1,dum2,dum3,ne]=randfilt(conv(nb,awb),conv(fsa,awa),0);
        end
        %% use the supplied noise signal
    else                                    % noise signal supplied
        hov=round(0.1*fs);                  % fading overlap width
        moff=min(2*any(m=='b')+any(m=='O'),2);      % moff=0='o',1='O', 2='b'
        mwr=min(2*any(m=='S')+any(m=='w'),2);       % mwr= 0='W',1='w', 2='S'
        nn=numel(nb);
        if size(nn,2)>1
            error('Noise signal must be a column vector');
        end
        if nargin<6
            fsa=fs;         % default sample rate is same as for speech
        end
        if fsa~=fs          % need to do resampling
            frat=fs/fsa;
            nbadd=ceil(10*max(frat,1)+1);         % half-length of filter at output sample rate
            nreq=ceil((ns+4*nbadd+hov+10)/frat);  % worst case number of input samples we need
            if nargout<3 && nreq<nn                              % just resample the bit we need
                ioff=0;                             % initial offset in original noise signal
                if moff<2
                    if mwr<2 && moff==0             % any offset allowed
                        ioff=floor(rand(1)*nn);
                    else
                        ioff=floor(rand(1)*(nn-nreq)); % choose offset to prevent wrapping
                    end
                end
                if nreq+ioff<nn % no wrapping required
                    nb=resample(reshape(nb(ioff+1:ioff+nreq),[],1),round(fs*10),round(fsa*10));
                    moff=2; % start at beginning
                    mwr=2;  % no wrapping
                else % resample portions at the beginning and end and define offset to avoid the middle
                    nb=resample([nb(1:nreq+ioff-nn+ceil(10/frat));nb(ioff+1:nn)],round(fs*10),round(fsa*10));
                    joff=length(nb)-floor((nn-ioff)*frat)+5-nbadd; % offset to use when wrapping
                    moff=3;
                end
            else % resample the entire noise signal
                nb=resample(nb(:),round(fs*10),round(fsa*10));
                if length(nb)<=2*nbadd
                    error('noise signal too short after resampling');
                end
            end
            nb=nb(nbadd+1:end-nbadd); % trim ends to avoid resampling transients
            nn=numel(nb);
        end
        %% determine the starting offset within the noise signal
        ioff=0;                         % initial offset in noise signal
        switch mwr
            case 2                      % mwr=2: no wrapping
                if nn<ns
                    error('noise signal too short');
                end
                if moff<2 % choose an offset
                    ioff=floor(rand(1)*(nn-ns));
                end
            case 1                      % mwr=1: abrupt wrapping
                switch moff
                    case 3
                        ioff=joff;  % use pre-calculated offset
                    case 1          % minimize number of wraps
                        ioff=floor(rand(1)*(nn-mod(ns-1,nn)));
                    case 0          % allow extra wraps
                        ioff=floor(rand(1)*nn);
                end

            case 0                   % mwr=0: faded wrapping
                hov=round(0.1*fs); % overlap is 0.1 seconds
                if nn<2*hov
                    error('noise signal too short');
                end
                switch moff
                    case 3
                        ioff=joff;  % use pre-calculated offset
                    case 1          % minimize number of wraps
                        ioff=floor(rand(1)*(nn+1-2*hov-mod(ns-hov-1,nn-hov)));
                    case 0          % allow any number of wraps
                        ioff=floor(rand(1)*nn-hov);
                end
                if nargin>=3 || ioff+ns>nn    	% we need to create a wrapped noise signal
                    wf=sin(0.25*pi*(1+cos((1:hov)*pi/(1+hov))))';
                    nb(nn-hov+1:nn)=nb(1:hov).*wf(hov:-1:1)+nb(nn-hov+1:nn).*wf;
                    nb(1:hov)=[];
                end
        end
        if any(m=='Z')
            ne=1;                       % noise power given as 0 dB
        elseif any(m=='A')
            sf=filter(awb,awa,nb);      % apply A-weighting
            if any(m=='P')
                ne=activlev(sf,fs);     % use P.56 for noise power
            else
                ne=mean(sf.^2);         % or just normal energy
            end
        else
            if any(m=='P')
                ne=activlev(nb,fs);     % use P.56 for noise power
            else
                ne=mean(nb.^2);          % or just normal energy
            end
        end
    end
    %% Determine scaling factors
    if any(m=='D')
        snre=snr;
    else
        snre=10^(0.1*snr); % convert from dB
    end
    if snre>1 % positive SNR -> fix the signal level
        if any(m=='n')
            sze=1;
        elseif any(m=='N')
            sze=snre;
        elseif any(m=='K')
            sze=ne*snre;
        elseif any(m=='k')
            sze=se;
        else
            sze=1/(1+1/snre);
        end
        nze=sze/snre;
    else % negative SNR -> fix the noise level
        if any(m=='n')
            nze=1/snre;
        elseif any(m=='N')
            nze=1;
        elseif any(m=='K')
            nze=ne;
        elseif any(m=='k')
            nze=se/snre;
        else
            nze=1/(1+snre);
        end
        sze=nze*snre;
    end
    pe=[se ne sze nze]; % powers
    if (se==0 && sze>0) || (ne==0 && nze>0) || sze==Inf || nze==Inf
        if (se==0 && sze>0) || sze==Inf
            error('Infinite gain for signal with mode ''%s'': input=%.1f dB, output=%.1f dB',m,db(se)/2,db(sze)/2);
        else
            error('Infinite gain for noise with mode ''%s'': input=%.1f dB, output=%.1f dB',m,db(ne)/2,db(nze)/2);
        end
    end
    gm=sqrt([sze/(se+(se==0)) nze/(ne+(ne==0))]); % gains for speech and noise
    p=pe;
    if ~any(m=='D') && nargout~=1
        mk=pe~=Inf & pe~=0;
        p(mk)=10*log10(pe(mk));
        p(pe==0)=-Inf;
    end
end
if gm(2)>0          % only generate noise if it has a non-zero gain
    nn=length(nb);  % caluculate the length of the noise signal or filter numerator
    if genf         % we need to generate new noise
        if nn==1 && numel(fsa)==1 && nb==1 && fsa==1
            n=randn(ns,1);
        else
            [n,zof]=filter(nb,fsa,randn(ns,1),zof);
        end
    else % use supplied noise signal
        n=nb(1+mod(ioff+nsc:ioff+ns-1+nsc,nn));
    end
    if any(m=='x')
        z=[gm(1)*s(:) gm(2)*n(:)];
    else
        z=gm(1)*s(:)+gm(2)*n(:);
    end
elseif any(m=='x')
    z=[gm(1)*s(:) zeros(ns,1)];
else
z=gm(1)*s(:);
end
%% create output state if required
if nargout>2
    fso.fs=fs;          % sample frequency
    fso.genf=genf;      % filtered white noise flag
    fso.nsc=nsc+ns; 	% cumulative sample count
    fso.nb=nb;            % noise signal or filter numerator
    fso.gm=gm; % gains for speech and noise
    fso.p=p; % output powers
    if genf
        fso.zof=zof;        % noise generating filter state
        fso.fsa=fsa;        % filter denominator
    else
        fso.ioff=ioff;      % offset in noise signal
    end
end
%% now plot if no output arguments
if ~isstruct(fsx) && (~nargout || any(m=='g'))
    tax=(1:ns)/fs;
    subplot(3,1,3)
    zp=sqrt(mean(z.^2));
    if zp>0
        lg=20*log10(zp);
    else
        lg=-Inf;
    end
    plot(tax,z,'-y',tax([1 ns]),[0 0],':k',tax([1 ns]),[zp zp],'-b');
    texthvc(tax(1),zp,sprintf(' %.1f dB',lg),'lbb');
    axisenlarge([-1.005 -1.05]);
    ylim=get(gca,'ylim');
    if snre==0
        lg=-Inf;
    elseif snre==Inf
        lg=Inf;
    else
        lg=10*log10(snre);
    end
    texthvc(tax(ns),ylim(2),sprintf('%.1f dB SNR ',lg),'rtb');
    xlabel('Time (s)');
    ylabel('Output');
    subplot(3,1,2)
    if nze==0
        plot(tax([1 ns]),[0 0],'-y');
        axisenlarge([-1.005 1]);
    else
        plot(tax,n,'-y',tax([1 ns]),[0 0],':k',tax([1 ns]),sqrt(pe(2))*[1 1],'-b');
        if pe(2)==0
            lg=-Inf;
        elseif pe(2)==Inf
            lg=Inf;
        else
            lg=10*log10(pe(2));
        end
        texthvc(tax(1),sqrt(pe(2)),sprintf(' %.1f dB',lg),'lbb');
        axisenlarge([-1.005 -1.05]);
        ylim=get(gca,'ylim');
        if pe(4)/pe(2)==0
            lg=-Inf;
        elseif pe(4)/pe(2)==Inf
            lg=Inf;
        else
            lg=10*log10(pe(4)/pe(2));
        end
        texthvc(tax(ns),ylim(2),sprintf('\\times %.1f dB ',lg),'rtb');
    end
    ylabel('Noise');
    subplot(3,1,1)
    plot(tax,s,'-y',tax([1 ns]),[0 0],':k',tax([1 ns]),sqrt(pe(1))*[1 1],'-b');
    if pe(1)==0
        lg=-Inf;
    elseif pe(1)==Inf
        lg=Inf;
    else
        lg=10*log10(pe(1));
    end
    texthvc(tax(1),sqrt(pe(1)),sprintf(' %.1f dB',lg),'lbb');
    axisenlarge([-1.005 -1.05]);
    ylim=get(gca,'ylim');
    if pe(3)/pe(1)==0
        lg=-Inf;
    elseif pe(3)/pe(1)==Inf
        lg=Inf;
    else
        lg=10*log10(pe(3)/pe(1));
    end
    texthvc(tax(ns),ylim(2),sprintf('\\times %.1f dB ',lg),'rtb');
    ylabel('Signal');
end
