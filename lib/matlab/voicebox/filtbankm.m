function [x,cf,il,ih]=filtbankm(p,n,fs,fl,fh,w)
%FILTBANKM determine matrix for a linear/mel/erb/bark-spaced filterbank [X,IL,IH]=(P,N,FS,FL,FH,W)
%
% Usage:
% (1) Calcuate the Mel-frequency Cepstral Coefficients
%
%        f=rfft(s);			              % rfft() returns only 1+floor(n/2) coefficients
%		 x=filtbankm(p,n,fs,0,fs/2,'m');  % n is the fft length, p is the number of filters wanted
%		 z=log(x*abs(f).^2);              % multiply x by the power spectrum
%		 c=dct(z);                        % take the DCT
%
% (2) Calcuate the Mel-frequency Cepstral Coefficients efficiently
%
%        f=fft(s);                                   % n is the fft length, p is the number of filters wanted
%        [x,cf,na,nb]=filtbankm(p,n,fs,0,fs/2,'m');  % na:nb gives the fft bins that are needed
%        z=log(x*(f(na:nb)).*conj(f(na:nb)));        % multiply x by the power spectrum
%		 c=dct(z);                                   % take the DCT
%
% (3) Plot the calculated filterbanks as a graph
%
%        plot((0:floor(n/2))*fs/n,filtbankm(p,n,fs,0,fs/2,'m')')   % fs=sample frequency
%
% (4) Plot the calculated filterbanks as a spectrogram
%
%        filtbankm(p,n,fs,0,fs/2,'m');
%
% Inputs:
%       p   number of filters in filterbank or the filter spacing in k-mel/bark/erb [ceil(4.6*log10(fs))]
%		n   length of fft
%           or [nfrq dfrq frq1] nfrq=number of input frequency bins, frequency increment (Hz), first bin freq (Hz)
%		fs  sample rate in Hz
%		fl  low end of the lowest filter in Hz (see 'h' option) [default = 0 or 30Hz for 'l' option]
%		fh  high end of highest filter in Hz [default = fs/2]
%		w   any sensible combination of the following:
%
%             'b' = bark scale instead of mel
%             'e' = erb-rate scale
%             'l' = log10 Hz frequency scale
%             'f' = linear frequency scale [default]
%             'm' = mel frequency scale
%
%             'n' - round to the nearest FFT bin so each row of x contains only one non-zero entry
%
%             'c' = fl & fh specify centre of low and high filters instead of edges
%             'h' = fl & fh are in mel/erb/bark/log10 instead of Hz
%             'H' = cf outputs are in mel/erb/bark/log10 instead of Hz
%
%		      'y' = lowest filter remains at 1 down to 0 frequency and
%			        highest filter remains at 1 up to nyquist freqency
%		            The total power in the fft is preserved (unless 'u' is
%		            specified).
%             'Y' = extend only at low frequency end (or high end if 'y' also specified)
%
%             'p' = input P specifies the number of filters [default if P>=1]
%             'P' = input P specifies the filter spacing [default if P<1]
%
%             'u' = input and output are power per Hz instead of power.
%             'U' = input is power but output is power per Hz.
%
%             's' = single-sided input: do not include symmetric negative frequencies (i.e. non-DC inputs have been doubled)
%             'S' = single-sided output: do not mirror the non-DC filter characteristics (i.e. double non-DC outputs)
%
%             'g' = plot filter coefficients as graph
%             'G' = plot filter coefficients as image [default if no output arguments present]
%
%
% Outputs:	x     a sparse matrix containing the filterbank amplitudes
%		          If the il and ih outputs are given then size(x)=[p,ih-il+1]
%                 otherwise size(x)=[p,1+floor(n/2)]
%                 Note that the peak filter values equal 2 to account for the power
%                 in the negative FFT frequencies.
%           cf    the filterbank centre frequencies in Hz (see 'H' option)
%		    il    the lowest fft bin with a non-zero coefficient
%		    ih    the highest fft bin with a non-zero coefficient
%
% The routine performs interpolation of the input spectrum by convolving the power spectrum
% with a triangular filter and then simulates a filterbank with asymetric triangular filters.

%
% References:
%
% [1] S. S. Stevens, J. Volkman, and E. B. Newman. A scale for the measurement
%     of the psychological magnitude of pitch. J. Acoust Soc Amer, 8: 185–19, 1937.
% [2] S. Davis and P. Mermelstein. Comparison of parametric representations for
%     monosyllabic word recognition in continuously spoken sentences.
%     IEEE Trans Acoustics Speech and Signal Processing, 28 (4): 357–366, Aug. 1980.

% Bugs/Suggestions
% (1) default frequencies won't work if the h option is specified
% (2) low default frequency is invalid if the 'l' option is specified

%      Copyright (C) Mike Brookes 1997-2009
%      Version: $Id: filtbankm.m 8483 2016-09-13 15:48:05Z dmb $
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

% Note "FFT bin_0" assumes DC = bin 0 whereas "FFT bin_1" means DC = bin 1

if nargin < 6
    w='f'; % default option: linear frequency scale
end
wr=' ';   % default warping is linear frequency
for i=1:length(w)
    if any(w(i)=='lebm');
        wr=w(i);
    end
end
if nargin < 5 || ~numel(fh)
    fh=0.5*fs; % max freq is the nyquist
end
if nargin < 4 || ~numel(fl)
    if wr=='l'
        fl=30;  % min freq is 30 Hz for log scale
    else
        fl=0; % min freq is DC
    end
end

fa=0;
if numel(n)>1
    nf=n(1);  % number of input frequency bins
    df=n(2);  % input frequency bin spacing
    if numel(n)>2
        fa=n(3); % frequency of first bin
    end
else
    nf=1+floor(n/2); % number of input frequency bins
    df=fs/n;  % input frequency bin spacing
end
fin0=fa+(0:nf-1)*df;  % input frequency bins

mflh=[fl fh];
if ~any(w=='h')             % convert Hz to mel/erb/...
    switch wr
        case 'm'
            mflh=frq2mel(mflh);       % convert frequency limits into mel
        case 'l'
            if fl<=0
                error('Low frequency limit must be >0 for l option');
            end
            mflh=log10(mflh);       % convert frequency limits into log10 Hz
        case 'e'
            mflh=frq2erb(mflh);       % convert frequency limits into erb-rate
        case 'b'
            mflh=frq2bark(mflh);       % convert frequency limits into bark
    end
end
melrng=mflh*(-1:2:1)';          % mel/erb/... range
% fn2=floor(n/2);     % bin index of highest positive frequency (Nyquist if n is even)
if isempty(p)
    p=ceil(4.6*log10(2*(fa+(nf-1)*df)));         % default number of output filters
end
puc=any(w=='P') || (p<1) && ~any(w=='p');
if any(w=='c')              % c option: specify fiter centres not edges
    if puc
        p=round(melrng/(p*1000))+1;
    end
    melinc=melrng/(p-1);
    mflh=mflh+(-1:2:1)*melinc;
else
    if puc
        p=round(melrng/(p*1000))-1;
    end
    melinc=melrng/(p+1);
end
%
% Calculate the FFT bins0 corresponding to the filters
%
cf=mflh(1)+(0:p+1)*melinc; % centre frequencies in mel/erb/... including dummy ends
cf(2:end)=max(cf(2:end),0); % only the first point can be negative
switch wr    % convert centre frequencies from mel/erb/... to Hz
    case 'l'
        mb=10.^(cf);
    case 'e'
        mb=erb2frq(cf);
    case 'b'
        mb=bark2frq(cf);
    case 'm'
        mb=mel2frq(cf);
    otherwise
        mb=cf;
end

% first sort out 2-sided input frequencies

fin=fin0;
fin(nf+1)=fin(nf)+df; % add on a dummy point at the high end
if fin(1)==0
    fin=[-fin(nf+1:-1:2) fin];
elseif fin(1)<=df/2
    fin=[-fin(nf+1:-1:1) fin];
elseif fin(1)<df
    fin=[-fin(nf+1:-1:1) fin(1)-df df-fin(1) fin];
elseif fin(1)==df
    fin=[-fin(nf+1:-1:1) 0 fin];
else
    fin=[-fin(nf+1:-1:1) df-fin(1) fin(1)-df fin];
end
nfin=length(fin);  % length of extended input frequency list

% now sort out the interleaving

fout=mb;  % output frequencies in Hz
lowex=any(w=='y')~=any(w=='Y');   % extend to 0 Hz
highex=any(w=='y') && (fout(end-1)<fin(end));  % extend at high end
if lowex
    fout=[0 0 fout(2:end)];
end
if highex
    fout=[fout(1:end-1) fin(end) fin(end)];
end
mfout=length(fout);
if any(w=='u') || any(w=='U')
    gout=fout(3:mfout)-fout(1:mfout-2);
    gout=2*(gout+(gout==0)).^(-1); % Gain of output triangles
else
    gout=ones(1,mfout-2);
end
if any(w=='S')
    msk=fout(2:mfout-1)~=0;
    gout(msk)=2*gout(msk); % double non-DC outputs for a 1-sided output spectrum
end
if any(w=='u')
    gin=ones(1,nfin-2);
else
    gin=fin(3:nfin)-fin(1:nfin-2);
    gin=2*(gin+(gin==0)).^(-1); % Gain of input triangles
end
msk=fin(2:end-1)==0;
if any(w=='s')
    gin(~msk)=0.5*gin(~msk); % halve non-DC inputs to change back to a 2-sided spectrum
end
if lowex
    gin(msk)=2*gin(msk);  % double DC input to preserve its power
end
foutin=[fout fin];
nfall=length(foutin);
wleft=[0 fout(2:mfout)-fout(1:mfout-1) 0 fin(2:nfin)-fin(1:nfin-1)]; % left width
wright=[wleft(2:end) 0]; % right width
ffact=[0 gout 0 0 gin(1:min(nf,nfin-nf-2)) zeros(1,max(nfin-2*nf-2,0)) gin(nfin-nf-1:nfin-2) 0]; % gain of triangle posts
% ffact(wleft+wright==0)=0; % disable null width triangles shouldn't need this if all frequencies are distinct
[fall,ifall]=sort(foutin);
jfall=zeros(1,nfall);
infall=1:nfall;
jfall(ifall)=infall; % unsort->sort index
ffact(ifall([1:max(jfall(1),jfall(mfout+1))-2 min(jfall(mfout),jfall(nfall))+2:nfall]))=0;  % zap nodes that are much too small/big

nxto=cumsum(ifall<=mfout);
nxti=cumsum(ifall>mfout);
nxtr=min(nxti+1+mfout,nfall);  % next input node to the right of each value (or nfall if none)
nxtr(ifall>mfout)=1+nxto(ifall>mfout); % next post to the right of opposite type (unsorted indexes)
nxtr=nxtr(jfall);  % next post to the right of opposite type (converted to unsorted indices) or if none: nfall/(mfout+1)

% each triangle is "attached" to the node at its extreme right end
% the general result for integrating the product of two trapesiums with
% heights (a,b) and (c,d) over a width x is (ad+bc+2bd+2ac)*w/6
%
% integrate product of lower triangles

msk0=(ffact>0);
msk=msk0 & (ffact(nxtr)>0); % select appropriate triangle pairs (unsorted indices)
ix1=infall(msk); % unsorted indices of leftmost post of pair
jx1=nxtr(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix1)-foutin(jx1-1); % length of right triangle to the left of the left post
yx=min(wleft(ix1),vfgx); % integration length
wx1=ffact(ix1).*ffact(jx1).*yx.*(wleft(ix1).*vfgx-yx.*(0.5*(wleft(ix1)+vfgx)-yx/3))./(wleft(ix1).*wleft(jx1)+(yx==0));

% integrate product of upper triangles

nxtu=max([nxtr(2:end)-1 0],1);
msk=msk0 & (ffact(nxtu)>0);
ix2=infall(msk); % unsorted indices of leftmost post of pair
jx2=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix2+1)-foutin(jx2); % length of left triangle to the right of the right post
yx=min(wright(ix2),vfgx); % integration length
yx(foutin(jx2+1)<foutin(ix2+1))=0; % zap invalid triangles
wx2=ffact(ix2).*ffact(jx2).*yx.^2.*((0.5*(wright(jx2)-vfgx)+yx/3))./(wright(ix2).*wright(jx2)+(yx==0));

% integrate lower triangle and upper triangle that ends to its right

nxtu=max(nxtr-1,1);
msk=msk0 & (ffact(nxtu)>0);
ix3=infall(msk); % unsorted indices of leftmost post of pair
jx3=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix3)-foutin(jx3); % length of upper triangle to the left of the lower post
yx=min(wleft(ix3),vfgx); % integration length
yx(foutin(jx3+1)<foutin(ix3))=0; % zap invalid triangles
wx3=ffact(ix3).*ffact(jx3).*yx.*(wleft(ix3).*(wright(jx3)-vfgx)+yx.*(0.5*(wleft(ix3)-wright(jx3)+vfgx)-yx/3))./(wleft(ix3).*wright(jx3)+(yx==0));

% integrate upper triangle and lower triangle that starts to its right

nxtu=[nxtr(2:end) 1];
msk=msk0 & (ffact(nxtu)>0);
ix4=infall(msk); % unsorted indices of leftmost post of pair
jx4=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix4+1)-foutin(jx4-1); % length of upper triangle to the left of the lower post
yx=min(wright(ix4),vfgx); % integration length
wx4=ffact(ix4).*ffact(jx4).*yx.^2.*(0.5*vfgx-yx/3)./(wright(ix4).*wleft(jx4)+(yx==0));

% now create the matrix

iox=sort([ix1 ix2 ix3 ix4;jx1 jx2 jx3 jx4]);
% [iox;[wx1 wx2 wx3 wx4]>0 ]
msk=iox(2,:)<=(nfall+mfout)/2;
iox(2,msk)=(nfall+mfout+1)-iox(2,msk);  % convert negative frequencies to positive
if highex
    iox(1,iox(1,:)==mfout-1)=mfout-2; % merge highest two output nodes
end
if lowex
    iox(1,iox(1,:)==2)=3; % merge lowest two output nodes
end

x=sparse(iox(1,:)-1-lowex,max(iox(2,:)-nfall+nf+1,1),[wx1 wx2 wx3 wx4],p,nf);
%
% sort out the output argument options
%
if ~any(w=='H')
    cf=mb;         % output Hz instead of mel/erb/...
end
cf=cf(2:p+1);  % remove dummy end frequencies
il=1;
ih=nf;
if any(w=='n') % round outputs to the centre of gravity bin
    sx2=sum(x,2);
    msk=full(sx2~=0);
    vxc=zeros(p,1);
    vxc(msk)=round((x(msk,:)*(1:nf)')./sx2(msk));
    x=sparse(1:p,vxc,sx2,p,nf);
end
if nargout > 2
    msk=full(any(x>0,1));
    il=find(msk,1);
    if ~numel(il)
        ih=1;
    elseif nargout >3
        ih=find(msk,1,'last');
    end
    x=x(:,il:ih);
end
if any(w=='u')
    sx=sum(x,2);
    x=x./repmat(sx+(sx==0),1,size(x,2));
end
%
% plot results if no output arguments or g option given
%
if ~nargout || any(w=='g') || any(w=='G') % plot idealized filters
    if ~any(w=='g') && ~any(w=='G')
        w=[w 'G'];
    end
    newfig=0;
    if  any(w=='g')
        plot(fa-df+(il:ih)*df,x');
        title(['filtbankm: mode = ' w]);
        xlabel(['Frequency (' xticksi 'Hz)']);
        ylabel('Weight');
        newfig=1;
    end
    
    if  any(w=='G')
        if newfig
            figure;
        end
        imagesc(fa-df+(il:ih)*df,1:p,x);
        axis 'xy'
        colorbar;
        cblabel('Weight');
        switch wr
            case 'l'
                type='Log-spaced';
            case 'e'
                type='Erb-spaced';
            case 'b'
                type='Bark-spaced';
            case 'm'
                type='Mel-spaced';
            otherwise
                type='Linear-spaced';
        end
        ylabel([type ' Filter']);
        xlabel(['Frequency (' xticksi 'Hz)']);
        title(['filtbankm: mode = ' w]);
    end
    
end