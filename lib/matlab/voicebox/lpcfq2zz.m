function zz=lpcfq2zz(f,q)
%LPCFQ2ZZ Convert frequencies and q factors to z-plane poles ZZ=(F,Q)
%all input values are in normalized Hz
% roots are at exp(2*pi*f*(-1/(2q) +- j)
% if f has more columns than q, remaining columns are real roots at -f

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: lpcfq2zz.m 713 2011-10-16 14:45:43Z dmb $
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

[nf,pf]=size(f);
if nargin < 2
   pq=0;
else
   pq=size(q,2);
end;
zz=zeros(nf,pf+pq);
if pq
   ii=1:pq;
   zz(:,2*ii-1)=exp(pi*f(:,ii).*(2i-q.^(-1)));
   zz(:,2*ii)=conj(zz(:,2*ii-1));
end
if pf>pq
   ii=1+pq:pf;
   zz(:,ii+pq)= exp(-2*pi*f(:,ii));
end
