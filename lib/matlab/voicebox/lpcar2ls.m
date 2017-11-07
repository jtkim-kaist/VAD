function ls=lpcar2ls(ar)
%LPCAR2LS convert ar polynomial to line spectrum pair frequencies LS=(AR)
% output vector elements will be in range 0 to 0.5
% the returned vector will be of length p

% This routine is nowhere near as efficient as it might be


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcar2ls.m 713 2011-10-16 14:45:43Z dmb $
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
p = p1-1;
p2 = fix(p/2);
d=0.5/pi;

if rem(p,2)		% odd order
  for k=1:nf
    aa=[ar(k,:) 0];
    r = aa + fliplr(aa);
    q = aa - fliplr(aa);
    fr = sort(angle(roots(r)));
    fq = [sort(angle(roots(deconv(q,[1 0 -1])))); 0];
    f = [fr(p2+2:p+1).' ; fq(p2+1:p).'];
    f(p+1) = [];
    ls(k,:) = d*f(:).';
  end
else
  for k=1:nf
    aa=[ar(k,:) 0];
    r = aa + fliplr(aa);
    q = aa - fliplr(aa);
    fr = sort(angle(roots(deconv(r,[1 1]))));
    fq = sort(angle(roots(deconv(q,[1 -1]))));
    f = [fr(p2+1:p).' ; fq(p2+1:p).'];
    ls(k,:) = d*f(:).';
  end
end
