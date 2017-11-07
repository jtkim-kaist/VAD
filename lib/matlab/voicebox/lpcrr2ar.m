function [ar,e]=lpcrr2ar(rr);
%LPCRR2AR convert autocorrelation coefs to ar coefs [AR,E]=(RR)
%E is the residual energy

% could test e each time and remove rows when it gets small


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcrr2ar.m 713 2011-10-16 14:45:43Z dmb $
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

[nf,p1]=size(rr);
p=p1-1;
ar=ones(nf,p1);
ar(:,2) = -rr(:,2)./rr(:,1);
e = rr(:,1).*(ar(:,2).^2-1);
for n = 2:p
   k = (rr(:,n+1)+sum(rr(:,n:-1:2).*ar(:,2:n),2)) ./ e;
   ar(:,2:n) = ar(:,2:n)+k(:,ones(1,n-1)).*ar(:,n:-1:2);
   ar(:,n+1) = k;
   e = e.*(1-k.^2);
end
e=-e;
