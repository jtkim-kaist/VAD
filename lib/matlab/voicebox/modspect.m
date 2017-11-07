function [c,qq,ff,tt,po]=modspect(s,fs,m,nf,nq,p)
%MODSPECT Calculate the modulation spectrum of a signal C=(S,FS,W,NC,P,N,INC,FL,FH)
%
%
% Simple use:
%
% Inputs:
%     s	  speech signal
%     fs  sample rate in Hz (default 11025)
%     m   mode string containing any sensible combination of the following letters; these
%         options override any specified in argument p.
%
%             g     apply AGC to the input speech
%
%           r/n/m   time domain windows: rectangular/hanning/hamming [def: hamming]
%           T/N/M   mel/mod-mel domain windows: triangular/hanning/hamming [def: triangular]
%
%            a/c    mel filters act in the amplitude or complex domain [def: power]
%            p/l    mel outputs are in power or log power [def: amplitude]
%             A     mod-mel filters act in the power domain [def: power]
%            P/L    mod-mel outputs are in amplitude or log power [def: amplitude]
%             k     include DC component of modulation spectrum [not implemented]
%
%             F     Take a DCT in the mel direction
%             Z     include 0'th order mel-cepstral coef
%
%             f     Take a DCT in the mod-mel direction
%             z     include 0'th order mod-mel cepstral coef
%
%             D     include d/df coefficients in * per mel (not valid with 'F')
%             d     include d/dt coefficients in * per second
%
%             u     give frequency axes in mels and mod-mels rather than in Hz
%
%             s     plot the result
%             S     plot intermediate results also
%
%     nf  four-element vector giving:
%           nf(1) = number of mel bins (or increment in k-mel) [0.1]
%           nf(2) = min filterbank frequency [40 Hz]
%           nf(3) = max filterbank frequency [10kHz]
%           nf(4) = number of DCT coefficientss excl 0'th [25].
%         Omitted values use defaults.
%     nq  four-element vector giving
%           nq(1) = number of mel-mod bins (or increment in 10 mod-mel units)
%           nq(2) = min filterbank frequency [50 Hz]
%           nq(3) = max filterbank frequency [400 Hz]
%           nq(4) = number of DCT coefficients excl 0'th [15]
%         Omitted values use defaults. Note that the mel-mod frequency scale
%         is like the mel scale  but reduced by a factor of 100.
%         I.e. mod-mel(f)=mel(100 f)/100; thus 10Hz corresponds to 10 mod-mels.
%     p   parameter structure containing some or all of the parameters listed
%         at the start of the program code below. p may be given as the
%         last input parameter even if some or all of the inputs m, nf and nq
%         have been omitted.
%

%
%
% Outputs:
%           c[*,nf,nt]  modulation spectrum where nt the number of time frames and
%                       nf is the normally number of mel frequency bins (but reduced by 1 if
%                       'F' is specified without 'Z').
%                       The middle index is for the output coefficients
%                       which are concatenated in the order [base, delta-p, delta-t].
%                       The number of base coefficients is normally the
%                       number of mon-mel filters but will be reduced by 1
%                       if 'f' is specified but 'z' is not.
%           qq          centre frequencies of modulation bins before any DCT (Hz or mel if 'u' set)
%           ff          centre frequencies of freqiency bins before any DCT (Hz or mel if 'u' set)
%           tt          time axis (sample s(i) is at time i/fs)
%           po          a copy of the parameter structure that was actually used.
%
% Algorithm Outline [parameters are in brackets]:
%
%   (1) Apply AGC to the speech to preserve a constant RMS level [tagc,tagt;'g']
%   (2) Divide into overlapping frames [tinc,tovl], window [twin;'rnm'], zero-pad [tpad,trnd] and FFT
%   (3) Apply mel filters in power or amplitude domain [fmin,fmax,fbin,fwin,fpow;'TNMac'=nf]
%   (4) Regarding each filterbank output as a time-varying signal power or amplitude [ppow;'pl'],
%       divide into overlapping frames [pinc,povl], window [pwin;'rnm'], zero-pad [ppad,prnd] and FFT
%   (5) Apply mod-mel filters in power or amplitude domain [mmin,mmax,mbin,mwin,mpow;'TNMA'=nq]
%   (6) Take optional log [qpow;'PL']
%   (7) Optionally apply a DCT in the mel [qdcp;'F'] and/or mod-mel [qdcq;'f'] directions preserving
%       the only some low quefrency coefficients [ dncp, dzcp;'Z'=nf, dncq, dzcq;'z'=nq]
%   (8) Calculate derivatives over mel-frequency [ddep,ddnp;'D'] and/or over time [ddet,ddnt;'d']
%
% Notes: All logs are limitied in dynamic range [logr], output units can be Hz or mel [unit;'u']
%        and optional plots can be generated [dplt;'pP'].
%

%
% [1] S. D. Ewert and T. Dau. Characterizing frequency slectivity for envelope fluctuations.
%     J. Acoust Soc Amer, 108 (3): 1181–1196, Sept. 2000.
% [2] M. L. Jepsen, S. D. Ewert, and T. Dau. A computational model of human auditory signal processing and perception.
%     J. Acoust Soc Amer, 124 (1): 422–438, July 2008.
% [3] G. Kim, Y. Lu, Y. Hu, and P. C. Loizou. An algorithm that improves speech intelligibility
%     in noise for normal-hearing listeners.
%     J. Acoust Soc Amer, 126 (3): 1486–1494, Sept. 2009. doi: 10.1121/1.3184603.
% [4] J. Tchorz and B. Kollmeier. Estimation of the signal-to-noise ratio with amplitude
%     modulation spectrograms. Speech Commun., 38: 1–17, 2002.
% [5] J. Tchorz and B. Kollmeier. SNR estimation based on amplitude modulation analysis with
%     applications to noise suppression.
%     IEEE Trans Speech Audio Processing, 11 (3): 184–192, May 2003. doi: 10.1109/TSA.2003.811542.

% bugs:
%  (3) should scale the derivatives so they get the correct units (e.g. power per second or per mel)
%  (5) should take care of the scaling so that FFT size does not affect
%  results
%  (4) sort out quefrency scale for graphs
%  (5) could add an option to draw an algorithm flow diagram

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: modspect.m 713 2011-10-16 14:45:43Z dmb $
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
persistent PP
if isempty(PP)

    % Algorithm constants: the first letter indicates the signal domain as follows:
    %     t=time, f=frequency, p=mel freq, m=modulation, q=mel-modulation, d=mod-cepstral

    PP.logr=120;                                % maximum dynamic range before taking log  (dB)
    PP.tagc=0;                                  % 1 = do agc
    PP.tagt=2;                                  % agc time constant
    PP.tinc=0.25e-3;                            % time domain frame increment
    PP.tovl=16;                                 % time domain overlap factor
    PP.tpad=30e-3;                              % time domain minimum FFT length (s). must be > 1/min spacing
    % of mel filters in Hz
    PP.trnd=1;                                  % 1 to round up time domain FFT to the next power of 2
    PP.twin='m';                                % time domain window shape: m=hamming, n=hann
    PP.fpow='p';                                % mel filters act on a=amplitude or p=power or c=complex FFT coefs
    PP.fmin=40;                                 % minimum mel filterbank frequency
    PP.fmax=10000;                              % maximum mel filterbank frequency (0 for Nyquist)
    PP.fwin='t';                                % mel filterbank filter shape: t=triangle, m=hamming, n=hann
    PP.fbin=0.1;                                 % number of mel filters or target spacing in k-mel
    PP.ppow='a';                                % modulation signal is a=amplitude or p=power or l=log-power
    PP.pinc=16e-3;                              % modulation frame increment
    PP.povl=2;                                  % modulation frame overlap factor
    PP.pwin='m';                                % mel domain time-window shape
    PP.ppad=200e-3;                              % mel domain minimum FFT length (s). must be > 1/min spacing
    % of mod-mel filters in Hz
    PP.prnd=1;                                  % 1 to round up mel domain FFT to the next power of 2
    PP.mpow='p';                                % transform mel domain a=amplitude or p=power
    PP.mwin='t';                                % mel-mod filterbank filter shape: t=triangle, m=hamming, n=hann
    PP.mbin=15;                                 % Number of mel-mod filters
    PP.mmin=50;                                 % minimum modulation frequency
    PP.mmax=400;                                % maximum modulation frequency (0 for Nyquist)
    PP.qpow='a';                                % Use mel-modulation domain a=amplitude or p=power or l=log-power
    PP.qdcq=0;                                  % 1 = Take DCT in modulation direction
    PP.qdcp=0;                                  % 1 = Take DCT in frequency direction
    PP.dzcq=0;                                  % 1 = include 0'th coefficient in modulation direction (always included in derivatives)
    PP.dzcp=0;                                  % 1 = include 0'th coefficient in frequency direction
    PP.dncp=30;                                 % number of frequency coefficient (excl 0'th)
    PP.dncq=15;                                 % number of mel-mod coefficient (excl 0'th)
    PP.ddet=0;                                  % 1 = include delta-time coefficients
    PP.ddnt=2;                                  % (length of time filter - 1)/2
    PP.ddep=0;                                  % 1 = include delta-freq coefficient
    PP.ddnp=1;                                  % (length of freq filter - 1)/2
    PP.unit=0;                                  % 1 = give frequency axes in mels and mod-mels
    PP.dplt=0;                                  % plotting bitfield: 1=new fig,2=link axes, 4=base coefs, 8=d/dp, 16=d/dt,
    % 32=waveform, 64=spectrogram, 128=mel-spectrogram, 256=mod-spectrogram, 512 frequency chart
    PP.cent=0.2;                                % centile to check for log plotting
    PP.cchk=0.2;                                % fraction to require above the centile to allow linear plot
end
po=PP;              % initialize the parameter structure to the default values
switch nargin
    case 0
        % list all fields
        nn=sort(fieldnames(PP));
        cnn=char(nn);
        fprintf('%d Voicebox parameters:\n',length(nn));

        for i=1:length(nn);
            if ischar(PP.(nn{i}))
                fmt='  %s = %s\n';
            else
                fmt='  %s = %g\n';
            end
            fprintf(fmt,cnn(i,:),PP.(nn{i}));
        end
        return
    case 1
        error('no sample frequency specified');
    case 2
        p=[];
    case 3
        if isstruct(m)
            p=m;
            m='';
        else
            p=[];
        end
        nf=[];
        nq=[];
    case 4
        if isstruct(nf)
            p=nf;
            nf=[];
        else
            p=[];
        end
        nq=[];
    case 5
        if isstruct(nq)
            p=nq;
            nq=[];
        else
            p=[];
        end
end
if isstruct(p)      % copy specified fields into po structure
    nn=fieldnames(p);
    for i=1:length(nn)
        if isfield(po,nn{i})
            po.(nn{i})=p.(nn{i});
        end
    end
end

% now sort out the mode flags

for i=1:length(m)
    switch m(i)
        case 'g'
            po.tagc=1;
        case {'r','n','m'}
            po.twin=m(i);
            po.pwin=m(i);
        case {'T','N','M'}
            po.fwin=lower(m(i));
            po.mwin=po.fwin;
        case {'a','c'}
            po.fpow=m(i);
        case {'p','l'}
            po.ppow=m(i);
        case 'A'
            po.mpow=lower(m(i));
        case {'P','L'}
            po.qpow=lower(m(i));
        case 'f'
            po.qdcq=1;
        case 'F'
            po.qdcp=1;
        case 'z'
            po.dzcq=1;
        case 'Z'
            po.dzcp=1;
        case 'd'
            po.ddet=1;
        case 'D'
            po.ddep=1;
        case 'u'
            po.unit=1;
        case 's'
            po.dplt=bitor(po.dplt,31+512);
        case 'S'
            po.dplt=bitor(po.dplt,1023);
    end
end

if ~nargout
    po.dplt=bitor(po.dplt,4);   % just plot base coefficients
end
if ~isempty(nf)
    po.dncp=nf;
end
if ~isempty(nq)
    po.dncq=nq;
end

lnf=length(nf);
if lnf>=1
    po.fbin=nf(1);
    if lnf>=2
        po.fmin=nf(2);
        if lnf>=3
            po.fmax=nf(3);
            if lnf>=4
                po.dncp=nf(4);
            end
        end
    end
end
lnq=length(nq);
if lnq>=1
    po.mbin=nq(1);
    if lnq>=2
        po.mmin=nq(2);
        if lnq>=3
            po.mmax=nq(3);
            if lnq>=4
                po.dncq=nq(4);
            end
        end
    end
end

% finally eliminate potential incomaptibilities

po.ddep=po.ddep & (po.qdcp==0);         % never differentiate along DCT axis
po.dzcq=po.dzcq | (po.qdcq==0);         % include all coefficients if no DCT
po.dzcp=po.dzcp | (po.qdcp==0);         % include all coefficients if no DCT

logth=10^(-po.logr/10);
ns=length(s);
axhandle=[];
%
% first apply AGC
%
if po.tagc
    tau=po.tagt*fs;      % time constant for power filtering
    ax2=[1 -exp(-1/tau)];
    bx2=sum(ax2);
    sm=mean(s(1:min(ns,round(tau))).^2);            % initialize filter using average of 1 tau
    if sm>0
        s=s./sqrt(filter(bx2,ax2,s.^2,-ax2(2)*sm));    % normalize to RMS=1
    end
end
if bitand(po.dplt,32)
    if bitand(po.dplt,1)
        figure;
    end
    plot((1:ns)/fs,s);
    axhandle(end+1)=gca;
    title('Input after AGC');
    xlabel('Time (s)');
end
%
% % calculate power spectra
%
ni=ceil(po.tinc*fs);    % frame increment in samples
nw=ni*po.tovl;           % window length
nf=ceil(max(nw,po.tpad*fs));          % FFT length: pad with zeros
if po.trnd
    nf=2^nextpow2(nf);  % round up FFT length to next power of 2
else
    nf=round(nf);
end
switch po.twin
    case 'm'
        w=hamming(nw+1)'; w(end)=[];        % Hamming window ([2] uses Hanning)
    case 'n'
        w=hanning(nw-1)';           % Hanning window: underlying period is nw
    case 'r'
        w=ones(nw,1);                 % rctangular window
end
fx=rfft(enframe(s,w,ni),nf,2);
[nft,nff]=size(fx);                 % number of frames and frequency bins
if bitand(po.dplt,64)
    if bitand(po.dplt,1)
        figure;
    end
    fp=fx.*conj(fx);
    imagesc(((0:nft-1)*ni+(nw+1)/2)/fs,(0:nff-1)*fs*0.001/nf,10*log10(max(fp,max(fp(:))*1e-6))');
    axhandle(end+1)=gca;
    hanc= colorbar;
    set(get(hanc,'ylabel'),'String','dB');
    axis('xy');
    title('Input Spectrogram');
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');
end
[mbk,ffm]=melbankm(po.fbin,nf,fs,po.fmin/fs,min(po.fmax/fs+(po.fmax==0),0.5),po.fwin);
switch [po.fpow po.ppow]
    case 'aa'
        px=abs(fx)*mbk';         % amplitude domain filtering
    case {'ap','al'}
        px=(abs(fx)*mbk').^2;         % amplitude domain filtering
    case 'pa'
        px=sqrt((fx.*conj(fx))*mbk');         % power domain filtering
    case {'pp','pl'}
        px=(fx.*conj(fx))*mbk';         % power domain filtering
    case 'ca'
        px=abs(fx*mbk');         % complex domain filtering
    case {'cp','cl'}
        px=fx*mbk';         % complex domain filtering
        px=px.*conj(px);   % convert to power
end
if po.ppow=='l'
    px=log(max((px),max(px(:))*logth)); % clip before log
end
if bitand(po.dplt,128)
    if bitand(po.dplt,1)
        figure;
    end
    switch po.ppow
        case 'a'
            imagesc(((0:nft-1)*ni+(nw+1)/2)/fs,ffm/1000,20*log10(max(px,max(px(:))*1e-3))');
        case 'p'
            imagesc(((0:nft-1)*ni+(nw+1)/2)/fs,ffm/1000,10*log10(max(px,max(px(:))*1e-6))');
        case 'l'
            imagesc(((0:nft-1)*ni+(nw+1)/2)/fs,ffm/1000,max(px*10/log(10),max(px(:))*10/log(10)-60)');
    end
    axhandle(end+1)=gca;
    hanc= colorbar;
    set(get(hanc,'ylabel'),'String','dB');
    axis('xy');
    title('Mel Spectrogram');
    xlabel('Time (s)');
    ylabel('Frequency (k-mel)');
end
npf=length(ffm);     % number of mel filters
mni=ceil(po.pinc*fs/ni);    % frame increment in spectral frames
mnw=mni*po.povl;        % window length
mnf=ceil(max(mnw,po.ppad*fs/ni));          % FFT length: pad with zeros
if po.prnd
    mnf=2^nextpow2(mnf);  % round up FFT length to next power of 2
else
    mnf=round(mnf);
end
nmt=fix((nft-mnw+mni)/mni); % number of modulation spectrum frames
ix=repmat((1:mnw)',1,nmt)+repmat((0:nmt-1)*mni,mnw,1);  % time index for all modumation spectrum frames
mx=rfft(reshape(px(ix(:),:),mnw,nmt*npf),mnf,1);          % find modulation spectrum
% nmm=size(mx,1);             % number of modulation spectrum bins [not needed]
[qbk,qqm]=melbankm(po.mbin,mnf,100*fs/ni,po.mmin*ni/fs,min(po.mmax*ni/fs+(po.mmax==0),0.5),po.mwin);    % multiply frq by 100 to make it alsmost log
nqq=length(qqm);
switch po.mpow
    case 'a'
        qx=reshape(qbk*abs(mx),nqq,nmt,npf);
    case 'p'
        qx=reshape(qbk*(mx.*conj(mx)),nqq,nmt,npf);
end
switch [po.mpow po.qpow]
    case 'aa'
        qx=reshape(qbk*abs(mx),nqq,nmt,npf);         % amplitude domain filtering
    case {'ap','al'}
        qx=reshape(qbk*abs(mx),nqq,nmt,npf).^2;         % amplitude domain filtering + power out
    case 'pa'
        qx=sqrt(reshape(qbk*(mx.*conj(mx)),nqq,nmt,npf));         % power domain filtering + amp out
    case {'pp','pl'}
        qx=reshape(qbk*(mx.*conj(mx)),nqq,nmt,npf);         % power domain filtering
end
if po.qpow=='l'
    qx=log(max((qx),max(qx(:))*logth)); % clip before log
end
tt=((1+mnw*ni+nw-ni)/2+(0:nmt-1)*mni*ni)/fs; % time axis of output frames

if bitand(po.dplt,256) && (~bitand(po.dplt,4) || po.qdcq>0 || po.qdcp>0)
    if bitand(po.dplt,1)
        figure;
    end
    switch po.qpow
        case 'a'
            dqx=20*log10(max(qx,max(qx(:))*1e-3));
        case 'p'
            dqx=10*log10(max(qx,max(qx(:))*1e-6));
        case 'l'
            dqx=max(qx*10/log(10),max(qx(:))*10/log(10)-60);
    end
    dqx(end+1,:,:)=max(dqx(:));  % insert a border
    dqx=reshape(permute(dqx,[2,1,3]),nmt,(nqq+1)*npf);
    ffq=ffm(1)+((0:(nqq+1)*npf)-(nqq-1)/2)*(ffm(2)-ffm(1))/(nqq+1); % mel frequencies
    imagesc(tt,ffq/1000,dqx');
    axhandle(end+1)=gca;
    hanc= colorbar;
    set(get(hanc,'ylabel'),'String','dB');
    axis('xy');
    title('Modulation Spectrogram');
    xlabel('Time (s)');
    ylabel('Frequency (k-mel)');
end
ndq=nqq;    % number of coefficients to use in the q direction
if po.qdcq
    dx=reshape(rdct(reshape(qx,nqq,nmt*npf)),nqq,nmt,npf);    % take dct in q direction
    ndq=min(ndq,po.dncq+1);     % "+1" to include the 0'th coefficient
else
    dx=qx;
end
ndf=npf;    % number of coefficients to use in the p direction
if po.qdcp
    dx=permute(reshape(rdct(reshape(permute(qx,[3,1,2]),npf,nqq*nmt)),npf,nqq,nmt),[2 3 1]);    % take dct in p direction
    ndf=min(ndf,po.dncp+1);   % "+1" to include the 0'th coefficient
elseif po.ddep                  % calculate the frequency derivative
    nv=po.ddnp;
    vv=(-nv:nv)*-3/(nv*(nv+1)*(2*nv+1)*(ffm(2)-ffm(1)));
    dxdp=filter(vv,1,dx,[],3);
    dxdp=dxdp(:,:,[repmat(2*nv+1,1,nv) 2*nv+1:npf repmat(npf,1,nv)]);   % replicate filter outputs at ends
end
if po.ddet                      % calculate the time derivative
    nv=po.ddnt;
    vv=(-nv:nv)*-3/(nv*(nv+1)*(2*nv+1)*(tt(2)-tt(1)));
    dxdt=filter(vv,1,dx,[],2);
    dxdt=dxdt(:,[repmat(2*nv+1,1,nv) 2*nv+1:nmt repmat(nmt,1,nv)],:);   % replicate filter outputs at ends
end
nqj=ndq-(po.dzcq==0);
nqk=nqj+ndq*(po.ddep+po.ddet);
npk=ndf-(po.dzcp==0);
c=zeros(nqk,npk,nmt);
c(1:nqj,:,:)=permute(dx(1+ndq-nqj:ndq,:,1+ndf-npk:ndf),[1,3,2]);        % base coefficients
if bitand(po.dplt,512)
    if bitand(po.dplt,1)
        figure;
    end
    ffx=repmat(mel2frq(ffm(:)),1,nqq)/1000;
    qqx=repmat(mel2frq(qqm(:)')/100,npf,1);
    plot(qqx,ffx,'xb');
    axis([[qqx(1) qqx(end)]*[1.1 -0.1; -0.1 1.1] [ffx(1) ffx(end)]*[1.1 -0.1; -0.1 1.1]]);
    title('Frequency bin centres');
    xlabel(sprintf('%.1f : %.1f : %.1f mod-mel = Modulation Frequency (Hz)',qqm(1)/100,(qqm(2)-qqm(1))/100,qqm(end)/100));
    ylabel(sprintf('%.0f : %.0f : %.0f mel = Frequency (kHz)',ffm(1),ffm(2)-ffm(1),ffm(end)));
end
if bitand(po.dplt,4)
    if bitand(po.dplt,1)
        figure;
    end
    dqx=dx(1+ndq-nqj:ndq,:,1+ndf-npk:ndf);
    dqxmx=max(abs(dqx(:)));
    dqxge=dqx>=0;
    dbfact=2-(po.qpow=='p');        % 2=amplitude, 1=power
    dqxa=max(abs(dqx),dqxmx*10^(-6/dbfact));  % clip at -60 dB
    dqxmn=min(dqxa(:));
    dblab='';
    if(mean(dqxa(:)>dqxmn+po.cent*(dqxmx-dqxmn)))<po.cchk
        dboff=abs(round(dbfact*10*log10(dqxmn)));
        dblab='{\pm}dB';
        if ~all(dqxge(:))       % not all positive
            dqxa=dqxa/dqxmn;    % force log to be always positive
            if dqxmn>1 && dboff~=0
                dblab=sprintf('{\\pm}dB - %d',dboff);
            else
                dblab=sprintf('{\\pm}dB + %d',dboff);
            end
        else
            dblab='dB';
        end
        dqx=dbfact*10*log10(dqxa).*(2*dqxge-1);
    end
    dqx(end+1,:,:)=max(dqx(:));  % insert a border
    dqx=reshape(permute(dqx,[2,1,3]),nmt,(nqj+1)*npk);
    ffq=ffm(1)+((0:(nqj+1)*npf)-(nqj-1)/2)*(ffm(2)-ffm(1))/(nqj+1); % mel frequencies
    imagesc(tt,ffq/1000,dqx');
    axhandle(end+1)=gca;
    hanc= colorbar;
    set(get(hanc,'ylabel'),'String',dblab);
    axis('xy');
    if po.qdcq
        title('Modulation spectrum DCT');
    else
        title('Modulation spectrum');
    end
    xlabel('Time (s)');
    if po.qdcp
        ylabel('Quefrency (k-mel^{-1})');
    else
        ylabel('Frequency (k-mel)');
    end
end

if po.ddep
    c(nqj+1:nqj+ndq,:,:)=permute(dxdp(1:ndq,:,1+ndf-npk:ndf),[1,3,2]);  % add on p delta
    if bitand(po.dplt,8)
        if bitand(po.dplt,1)
            figure;
        end
        dqx=dxdp(1:ndq,:,1+ndf-npk:ndf);
        dqxmx=max(abs(dqx(:)));
        dqxge=dqx>=0;
        dbfact=2-(po.qpow=='p');        % 2=amplitude, 1=power
        dqxa=max(abs(dqx),dqxmx*10^(-6/dbfact));  % clip at -60 dB
        dqxmn=min(dqxa(:));
        dblab='';
        if(mean(dqxa(:)>dqxmn+po.cent*(dqxmx-dqxmn)))<po.cchk
            dboff=abs(round(dbfact*10*log10(dqxmn)));
            dblab='{\pm}dB';
            if ~all(dqxge(:))       % not all positive
                dqxa=dqxa/dqxmn;    % force log to be always positive
                if dqxmn>1 && dboff~=0
                    dblab=sprintf('{\\pm}dB - %d',dboff);
                else
                    dblab=sprintf('{\\pm}dB + %d',dboff);
                end
            else
                dblab='dB';
            end
            dqx=dbfact*10*log10(dqxa).*(2*dqxge-1);
        end
        dqx(end+1,:,:)=max(dqx(:));  % insert a border
        dqx=reshape(permute(dqx,[2,1,3]),nmt,(ndq+1)*npk);
        ffq=ffm(1)+((0:(ndq+1)*npf)-(ndq-1)/2)*(ffm(2)-ffm(1))/(ndq+1); % mel frequencies
        imagesc(tt,ffq/1000,dqx');
        axhandle(end+1)=gca;
        hanc= colorbar;
        set(get(hanc,'ylabel'),'String',dblab);
        axis('xy');
        if po.qdcq
            title('Modulation spectrum DCT = freq derivative');
        else
            title('Modulation spectrum = freq derivative');
        end
        xlabel('Time (s)');
        if po.qdcp
            ylabel('Quefrency (k-mel^{-1})');
        else
            ylabel('Frequency (k-mel)');
        end
    end
end
if po.ddet
    c(nqk-ndq+1:nqk,:,:)=permute(dxdt(1:ndq,:,1+ndf-npk:ndf),[1,3,2]);  % add on t delta
    if bitand(po.dplt,16)
        if bitand(po.dplt,1)
            figure;
        end
        dqx=dxdt(1:ndq,:,1+ndf-npk:ndf);
        dqxmx=max(abs(dqx(:)));
        dqxge=dqx>=0;
        dbfact=2-(po.qpow=='p');        % 2=amplitude, 1=power
        dqxa=max(abs(dqx),dqxmx*10^(-6/dbfact));  % clip at -60 dB
        dqxmn=min(dqxa(:));
        dblab='';
        if(mean(dqxa(:)>dqxmn+po.cent*(dqxmx-dqxmn)))<po.cchk
            dboff=abs(round(dbfact*10*log10(dqxmn)));
            dblab='{\pm}dB';
            if ~all(dqxge(:))       % not all positive
                dqxa=dqxa/dqxmn;    % force log to be always positive
                if dqxmn>1 && dboff~=0
                    dblab=sprintf('{\\pm}dB - %d',dboff);
                else
                    dblab=sprintf('{\\pm}dB + %d',dboff);
                end
            else
                dblab='dB';
            end
            dqx=dbfact*10*log10(dqxa).*(2*dqxge-1);
        end
        dqx(end+1,:,:)=max(dqx(:));  % insert a border
        dqx=reshape(permute(dqx,[2,1,3]),nmt,(ndq+1)*npk);
        ffq=ffm(1)+((0:(nqj+1)*npf)-(nqj-1)/2)*(ffm(2)-ffm(1))/(nqj+1); % mel frequencies
        imagesc(tt,ffq/1000,dqx');
        axhandle(end+1)=gca;
        hanc= colorbar;
        set(get(hanc,'ylabel'),'String',dblab);
        axis('xy');
        if po.qdcq
            title('Modulation spectrum DCT = time derivative');
        else
            title('Modulation spectrum = time derivative');
        end
        xlabel('Time (s)');
        if po.qdcp
            ylabel('Quefrency (k-mel^{-1})');
        else
            ylabel('Frequency (k-mel)');
        end
    end
end
if po.unit                      % frequency units are in mels
    ff=ffm;
    qq=qqm/100;
else
    ff=mel2frq(ffm);
    qq=mel2frq(qqm)/100;
end
if length(axhandle)>1
    if ~bitand(po.dplt,32+64)
        linkaxes(axhandle)
    else
        linkaxes(axhandle,'x')
    end
end


