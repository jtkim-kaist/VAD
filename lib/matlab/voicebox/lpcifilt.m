function u=lpcifilt(s,ar,t,dc,fade)
%LPCIFILT Apply inverse filter to speech signal U=(S,AR,T,DC,FADE)
%
% Inputs:	S             speech signal
%				AR(nf,p+1)    array of ar coefficients; one row per frame
%				T             column vector giving the index of the first sample in each frame
%           DC(nf)        optional dc component will be subtracted from the signal
%           FADE          AR coefficients will be linearly interpolated for FADE samples
%                         either side of frame boundaries
%

% Example usage: generate an inverse filtered waveform
%
%      [sp,fs]=readwav('infile.wav');
%      lpcord=2+floor(fs/1000);
%      spp=filter([1 -exp(-2*pi*50/fs)],1,sp);                    % preemphasis zero is at 50 Hz
%      [lpar,lpe,lpk]=lpcauto(spp,lpcord,floor([0.01 0.02]*fs));  % 20ms frame with 10ms frame increment
%      overlap=lpk(1,2)-lpk(2,1)+1;                               % overlap between adjacent frames
%      u=lpcifilt(sp,lpar,lpk(:,1)+floor(overlap/2),0,overlap/4); % do inverse filtering
%      writewav(u,fs,'outfile.wav');

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcifilt.m 2647 2013-01-23 11:05:41Z dmb $
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

[nf,p1]=size(ar);
if nargin<4 | isempty(dc)
   dc=zeros(nf,1);
elseif length(dc)==1
   dc=dc(ones(nf,1));
end
if nf==1
   u=filter(ar,1,s-dc);
else
   p=p1-1;
   ns=length(s);
   if nargin<5 | isempty(fade)
      fade=0;
      if nargin<3 | isempty(t)
         t=p1+(0:nf-1)*(ns-p)/nf;
      end
   end
   u=zeros(ns,1);
   if fade<1
      x0=(max(1,ceil(t(nf)-p)):ns)';
      u(x0)=filter(ar(nf,:),1,s(x0)-dc(nf));
      for i=nf-1:-1:2
         x0=(max(1,ceil(t(i)-p)):ceil(t(i+1)-1))';
         u(x0)=filter(ar(i,:),1,s(x0)-dc(i));
      end
      x0=(1:ceil(t(2)-1))';
      u(x0)=filter(ar(1,:),1,s(x0)-dc(1));
      
   else
      tb=min(t(nf)+fade,(t(nf)+ns)/2);
      ta=max(t(nf)-fade,(t(nf-1)+t(nf))/2);
      t0=max(1,ceil(ta)-p);
      x0=(t0:ns)';
      xb=x0-tb;
      u(x0)=filter(ar(nf,:),1,s(x0)-dc(nf)).*max((1+(xb-abs(xb))/(2*(tb-ta))),0);
      for i=nf-1:-1:2
         td=tb;
         tc=ta;
         tb=min(t(i)+fade,ta);
         ta=max(t(i)-fade,(t(i-1)+t(i))/2);
         t1=floor(td);
         t0=max(1,ceil(ta)-p);
         x0=(t0:t1)';
         xb=x0-tb;
         xc=x0-tc;
         u(x0)=u(x0)+filter(ar(i,:),1,s(x0)-dc(i)).*max((1+(xb-abs(xb))/(2*(tb-ta))-(xc+abs(xc))/(2*(td-tc))),0);
      end
      t1=floor(tb);
      x0=(1:t1)';
      xc=x0-ta;
      u(x0)=u(x0)+filter(ar(1,:),1,s(x0)-dc(1)).*(1-(xc+abs(xc))/(2*(tb-ta)));
      
   end
end

