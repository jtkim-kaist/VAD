function x=windinfo(w,fs)
%WINDINFO window information and figures of merit X=(W,FS)
%
%  Inputs:  W        is a vector containing the window
%           FS       is the sampling frequency (default=1)
%
% Outputs:  X.len         length of the window (samples)
%           X.nw          length of the window (samples)
%           X.ewgdelay    energy centroid delay from first sample (samples)
%           X.dcgain      DC gain (dB)
%           X.sidelobe    maximum sdelobe level in dB relative to DC gain
%           X.falloff     rate at which sidelobes decay (dB/octave)
%           X.enbw        equivalent noise bandwidth (*fs/len Hz)
%           X.scallop     scalloping loss (dB)
%           X.ploss       processing loss (dB)
%           X.wcploss     worst case processing loss (dB)
%           X.band3       3dB bandwidth (Hz)
%           X.band6       6 dB bandwidth (Hz)
%           X.band0       essential bandwidth (to first minimum) (Hz)
%           X.gain0       gain at first minimum (Hz)
%           X.olc50       50% overlap correction
%           X.olc75       75% overlap correction
%           X.cola        overlap factors giving constant overlap add (exlcuding multiples)
%           X.cola2       as X.cola but for squared window
%
% If no output argument is given, the window and frequency response
% will be plotted e.g. windinfo(windows('hamming',256,'ds'),256);
%
% To obtain the figures of merit listed in Table 1 of [1] set
% fs = length(W), multiply X.olc50 and X.olc75 by 100%. The "coherent gain
% listed in the table is 10^(x.dcgain/20)/(max(w)*length(w)).
%
%  [1]  F. J. Harris. On the use of windows for harmonic analysis with the
%       discrete fourier transform. Proc IEEE, 66 (1): 51–83, Jan. 1978.

%	   Copyright (C) Mike Brookes 2009-2014
%      Version: $Id: v_windinfo.m 6801 2015-09-12 09:30:42Z dmb $
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
%

if nargin<2
    fs=1;
end
w=w(:);
nw=length(w);
x.len=nw/fs;
x.nw=nw;
% energy weighted group delay = centre of energy
x.ewgdelay=((1:nw)*w.^2/sum(w.^2)-1)/fs;
% now calculate spectrum
of=16;      % spectrum oversample factor must be even
nwo=of*nw;
f=rfft(w,nwo);
p=f.*conj(f);
% sidelobe attenuation is maximum peak (note DC peak at p(1) is not found)
[kp,vp]=v_findpeaks(p,'q');
[kt,vt]=v_findpeaks(-p,'q');
if ~numel(kp)
    x.sidelobe=10*log10(min(p)/p(1));
else
    x.sidelobe=10*log10(max(vp)/p(1));
end
np=length(kp);
ipa=floor(np/4);
if ~ipa
    x.falloff=0;
else
    ipb=floor(np/2);
    x.falloff=10*log10(vp(ipb)/vp(ipa))/log2((ipb-1)/(ipa-1));
end
sumw2=sum(w.^2);
sumw=sum(w);
enbwbin=nw*sumw2/sumw^2;
x.enbw=enbwbin*fs/nw;
x.dcgain=20*log10(sumw);
% do linear interpolation in p() to find 3dB and 6dB points
p3=0.5*p(1);
i3=find(p<p3,1);
if ~numel(i3)
    x.band3=Inf;
    x.band6=Inf;
else
    x.band3=2*(i3-(p3-p(i3))/(p(i3-1)-p(i3))-1)/of*fs/nw;
    p6=0.25*p(1);
    i6=find(p<p6,1);
    x.band6=2*(i6-(p6-p(i6))/(p(i6-1)-p(i6))-1)/of*fs/nw;
end
% do linear interpolation in f() to find closest approach to the origin
if~numel(kt)
    x.band0=Inf;
    x.gain0=0;
else
    i0=floor(kt(1));
    df=f(i0+1)-f(i0);
    j0=-real(f(i0)*conj(df))/abs(df)^2;
    x.band0=2*(i0+j0-1)/of*fs/nw;
    p0=abs(f(i0)+j0*df)^2;
    if p0>0
        x.gain0=10*log10(p0/p(1));
    else
        x.gain0=-Inf;
    end
end
% overlap factors
i50=round(nw*0.5);
x.olc50=sum(w(1:nw-i50).*w(1+i50:nw))/sumw2;
i75=round(nw*0.25);
x.olc75=sum(w(1:nw-i75).*w(1+i75:nw))/sumw2;
% processing loss and scalloping loss
x.scallop=10*log10(p(1)/p(1+of/2));
x.ploss=10*log10(enbwbin);
x.wcploss=x.ploss+x.scallop;
co=zeros(1,nw);
co2=co;
w2=w.^2;
for i=1:nw
    if i*fix(nw/i)==nw  % check i is a factor of the window length
        if co(i)==0
            co(i)=all(abs(sum(reshape(w',nw/i,i),2)-i*mean(w))<max(abs(w))*1e-4);
            if co(i)
                co(2*i:i:nw)=-1; % ignore multiples
            end
        end
        if co2(i)==0
            co2(i)=all(abs(sum(reshape(w2',nw/i,i),2)-i*mean(w2))<max(w2)*1e-4);
            if co2(i)
                co2(2*i:i:nw)=-1; % ignore multiples
            end
        end
    end
end
co(co<0)=0;
co2(co2<0)=0;
x.cola=find(co);
x.cola2=find(co2);
%
% now plot it if no output arguments given
%
if ~nargout
    clf;

    subplot(212);
    nf=min(max(floor(2*max(x.band6,x.band0)*of*nw/fs)+1,of*8),length(p));
    ff=(0:nf-1)*fs/(of*nw);
    fqi=[x.enbw x.band3 x.band6]/2;
    if ff(end)>2000
        ff=ff/1000;
        fqi=fqi/1000;
        xlab='kHz';
    else
        xlab='Hz';
    end
    dbrange=min(100,-1.5*x.sidelobe);
    dd=10*log10(max(p(1:nf),p(1)*0.1^(dbrange/10)));
    ffs=[0 ff(end)];
    dbs=repmat(x.dcgain+x.sidelobe,1,2);
    ffb=[0 fqi(1) fqi(1)];
    dbb=[dd(1) dd(1) dd(1)-dbrange];
    ff3=[0 fqi(2) fqi(2)];
    db3=[dd(1)+db(0.5)/2 dd(1)+db(0.5)/2 dd(1)-dbrange];
    ff6=[0 fqi(3) fqi(3)];
    db6=[dd(1)+db(0.5) dd(1)+db(0.5) dd(1)-dbrange];
    area(ffb,dbb,max(dd)-dbrange,'facecolor',[1 0.7 0.7]);
    hold on
    plot(ffs,dbs,':k',ff3,db3,':k',ff6,db6,':k',ffb,dbb,'r',ff,dd,'b');
    legend(['Equiv Noise BW = ' sprintsi(x.enbw,-2) 'Hz'],['Max sidelobe = ' sprintf('%.0f',x.sidelobe) ' dB'],['-3 & -6dB BW = ' sprintsi(x.band3,-2) 'Hz & ' sprintsi(x.band6,-2) 'Hz']);
    hold off
    axis([0 ff(end) max(dd)-dbrange max(dd)+2]);
    ylabel('Gain (dB)');
    xlabel(sprintf('Freq (%s)',xlab));
    %
    % Now plot the window itself
    %
    subplot(211);
    tax=(0:nw-1)/fs-x.ewgdelay;
    area(tax,w);
    ylabel('Window');
    xlabel('Time (s)');
    dtax=(tax(end)-tax(1))*0.02;
    axv=[tax(1)-dtax tax(end)+dtax min(0,min(w)) max(w)*1.05];
    texthvc(tax(end),max(w),sprintf('N=%d',nw),'rtk');
    if length(x.cola)>3
        tcola=sprintf(',%d',x.cola(1:3));
        tcola=[tcola ',…'];
    else
        tcola=sprintf(',%d',x.cola);
    end
    if length(x.cola2)>3
        tcola2=sprintf(',%d',x.cola2(1:3));
        tcola2=[tcola ',…'];
    else
        tcola2=sprintf(',%d',x.cola2);
    end
    texthvc(tax(1),max(w),sprintf('COLA=%s\nCOLA^2=%s',tcola(2:end),tcola2(2:end)),'ltk');
    axis(axv);
end



