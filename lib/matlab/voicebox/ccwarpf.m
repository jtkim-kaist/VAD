function m=ccwarpf(f,n,s)
%CCWARPF  Warp cepstral coefficients M=(F,N,S) 
% f(1) is the original sample freq, f(2) is the new sample freq
% n(1) is the original number of coefficients, n(2) is the new number
% s is a string: s(1),s(2) =l for linear, m for mel frequency, use capitals if c0 included


%      Copyright (C) Mike Brookes 1998
%      Version: $Id: ccwarpf.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<3
   s='ll';
end
if length(f)<2
   f(2)=1;
end
if length(n)<2
   n(2)=n(1);
end
z=s<'a';
s=s+32*z;
if all(s=='l')
   k=1:n(2)-z(2);
   ff=((1:n(1)).'-z(1))*f(2)/f(1);
   fa=2*sin(ff*pi).*ff/pi;
   fb=ff.^2;
   ka=1-2*rem(k,2);
   kb=k.^2;
   r1=ones(n(1),1);
   c1=ones(1,n(2)-z(2));
   a=fa(:,c1).*ka(r1,:);
   b=fb(:,c1)-kb(r1,:);
   f0=find(fix(ff)==ff);
   if length(f0)
      a(f0,:)=ff(f0,c1)==k(ones(length(f0),1),:);
      b(f0,:)=1;
   end
   m=a./b;
   if z(2)
      m=[[1; 0.5*fa(2:n(1))./fb(2:n(1))] m];
   end
end



