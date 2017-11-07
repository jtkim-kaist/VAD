function s=lpcconv(from,to,x,y,np)
%LPCCONV(from,to,x,y)->s convert between LPC parameter sets
%
% The output is a string that may be passed to eval(s)
% x and y are optionally the input and output matrices
% and np the new value of the order p.
% with one frame stored per row. from and to are taken
% from the following list which also gives the column dimension:
%
%  1 ar p+1  Autoregressive coevfficients: ar(1)=1 always.
%  2 cc  p   Complex cepstral coefficients
%  3 ls  p   Line spectrum pair frequencies (normalized Hz)
%  4 zz  p   Z-plane roots
%  5 ss  p   S-plane roots (normalized Hz)
%  6 rf  p   Reflection coefficients (= -PARCOR coefs)
%  7 ao  p   Area ratios
%  8 aa p+2  Vocal tract areas: aa(p+2)=1 always
%  9 rr p+1  Autocorrelation coefficients
% 10 dl  p   DCT of log area function
% 11 lo  p   Log area ratios
% 12 la p+1  Log areas: la(1)=0 always
% 13 ra p+1  Autocorrelation coefs of inverse filter
% 14 ff p+2  Fourier transform of forward filter (all-pole)
% 15 pf p+2  Power spectrum of forward filter (all-pole)
% 16 gc  p   Gain and cos(w) of each formant
% 17 im p+1  Impulse response of forward filter

%	   Copyright (C) Mike Brookes 1998
%      Version: $Id: lpcconv.m 713 2011-10-16 14:45:43Z dmb $
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

nm=['aa';'am';'ao';'ar';'cc';'db';'dl';'ff';'fq';'im';'is';'la';'lo';'ls';'pf';'ra';'rf';'rr';'ss';'zz';];
nx=[...
      0 17  3 17 17 17  7 17  0 17 17 17 17 17 17 17 17 17 17 17;...
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
     17 17  0 17 17 17 17 17  0 17 17 17 17 17 17 17 17 17 17 17;...
     17 18 17  0  5  6 17  8  0 10 17 17 17 14 15 16 17 18 20 20;...
      4  4  4  4  5  6  4  4  0  4  4  4  4  4 15  4  4  4  4  4;...
     15 15 15 15 15  0 15 15  0 15 15 15 15 15 15 15 15 15 15 15;...
      1  1  1  1  1  1  0  1  0  1  1  1  1  1  1  1  1  1  1  1;...
     15 15 15 15 15 15 15  0  0 15 15 15 15 15 15 15 15 15 15 15;...
     20 20 20 20 20 20 20 20  0 20 20 20 20 20 20 20 20 20 20 20;...
      4  4  4  4  4  4  4  4  0  0  4  4  4  4  4  4  4  4  4  4;...
     17 17 17 17 17 17 17 17  0 17  0 17 17 17 17 17 17 17 17 17;...
     17 17 17 17 17 17 17 17  0 17 17  0 17 17 17 17 17 17 17 17;...
     17 17 17 17 17 17 17 17  0 17 17 17  0 17 17 17 17 17 17 17;...
      4  4  4  4  4  4  4  4  0  4  4  4  4  0  4  4  4  4  4  4;...
      5 18  5  5  5  5  5  5  0  5  5  5  5  5  0  5  5 18  5  5;...
     15 15 15 15 15 15 15 15  0 15 15 15 15 15 15  0 15 15 15 15;...
      1 18  3  4  4  4  1  4  0  4 11 12 13  4  4  4  0 18  4  4;...
      4  2  4  4  4  4  4  4  0  4  4  4  4  4  4  4  4  0  4  4;...
     20 20 20 20 20 20 20 20  0 20 20 20 20 20 20 20 20 20  0 20;...
      4  4  4  4  5  4  4  4  0  4  4  4  4  4  4  4  4  4 19  0;...
];
na=size(nm,1);
b=256*nm(:,1)+nm(:,2);
jf=find(b==256*from(1)+from(2));
jt=find(b==256*to(1)+to(2));
if length([jf jt])~=2
   [x,idx]=sort(b);
  error(sprintf('lpcxx2yy types are: %s',[nm(idx,:)';' '*ones(1,na)]));
end
if nargin<3 s=nm(jf,:); else s=x; end
while jf ~= jt
  jn=nx(jf,jt);
  if jn==0
     error(sprintf('cannot convert between %s and %s',nm(jf,:),nm(jt,:)));
 end
  s=sprintf('lpc%s2%s(%s)',nm(jf,:),nm(jn,:),s);
  jf=jn;
  end
if nargin<4 sn=nm(jt,:); else sn=y; end
s=sprintf('%s=%s;',sn,s);
    


