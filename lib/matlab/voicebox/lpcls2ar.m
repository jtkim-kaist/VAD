function ar=lpcls2ar(ls)
%LPCLS2AR convert line spectrum pair frequencies to ar polynomial AR=(LS)
% input vector elements should be in the range 0 to 0.5


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcls2ar.m 713 2011-10-16 14:45:43Z dmb $
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

[nf,p]=size(ls);
p1=p+1;
p2 = p1*2;
ar=zeros(nf,p1);
for k=1:nf
  le=exp(ls(k,:)*pi*2i);
  lf=[1 le -1 conj(fliplr(le))];
  y=real(poly(lf(1:2:p2)));
  x=real(poly(lf(2:2:p2)));
  ar(k,:)=(x(1:p1)+y(1:p1))/2;
end
