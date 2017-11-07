function pp=lpcar2pp(ar)
%LPCAR2PP LPC: Convert ar filter to power spectrum polynomial in cos(w) PP=(AR)
% the pp is a polynomial such that |polyval(ar,e^jw)|^2 = polyval(pp,cos(w))

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: lpcar2pp.m 713 2011-10-16 14:45:43Z dmb $
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

persistent tp;
[nf,p1]=size(ar);
ra=zeros(nf,p1);
for i=1:p1
   ra(:,i)=sum(ar(:,1:p1+1-i).*ar(:,i:p1),2);
end
% we check here if p is the same as before and reuse the tp matrix
if size(tp,1)~=p1
   p=p1-1;
   % chebyshev polynomials up to order p
   tp=zeros(p1,p1);
   tp(1,p1)=2;
   tp(2,p)=2;
   for i=3:p1
      tp(i,p+2-i:p)=2*tp(i-1,p+3-i:p1)-tp(i-2,p+2-i:p);
      tp(i,p1)=-tp(i-2,p1);
   end
   tp(1,p1)=1;
end
pp=ra*tp;
