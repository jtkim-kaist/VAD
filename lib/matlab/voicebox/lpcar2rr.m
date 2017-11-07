function rr=lpcar2rr(ar,p)
%LPCAR2RR Convert autoregressive coefficients to autocorrelation coefficients RR=(AR,P)
% The routine calculated the autocorrelation coefficients of the signal
% that results from feeding unit-variance, zero-mean noise into the all-pole filter
% Input:   ar(:,n+1)  Autoregressive coefficients including 0'th coefficient
% Output:  rr(:,p+1)    Autocorrelation coefficients including 0'th order coefficient
% If p is not specified it is taken to be n


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcar2rr.m 713 2011-10-16 14:45:43Z dmb $
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

k=ar(:,1).^(-2);
if size(ar,2)==1
   rr=k;
else
   if nargin>1
      rr=lpcrf2rr(lpcar2rf(ar),p).*k(:,ones(1,p+1));
   else
      rr=lpcrf2rr(lpcar2rf(ar)).*k(:,ones(1,size(ar,2)));
   end
end
