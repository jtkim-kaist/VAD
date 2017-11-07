function rr=lpcpf2rr(pf,p)
%LPCPF2RR convert power spectrum to autocorrelation coefs RR=(PF,P)
% Note that these will only be accurate if the power spectrum is much longer than p


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcpf2rr.m 2460 2012-10-29 22:20:45Z dmb $
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

[nf,p2]=size(pf);
if nargin<2 p=p2-2; end;
ir=irfft(pf,[],2);
if p>p2-2
  rr=[ir(:,1:p2-1) zeros(nf,p+2-p2)];
else
  rr=ir(:,1:p+1);
end
