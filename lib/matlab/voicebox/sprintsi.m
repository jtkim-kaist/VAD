function s=sprintsi(x,d,w)
%SPRINTSI Print X with SI multiplier S=(X,D,W)
% D is number of decimal places (+ve) or significant digits (-ve) [default=-3]
% |W| is total width including multiplier
% if W<=0 then trailing 0's will be eliminated
%
% Example: sprintsi(2345,-2) gives '2.3 k'

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: sprintsi.m 4966 2014-08-05 18:20:41Z dmb $
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

if nargin<3 w=0; end;
if nargin<2 d=-3; end;
f='afpnum kMGT';
e=max(-18,min(12,floor(log10(abs(x)))));
k=floor(e/3);
dp=max([0 d 3*k-d-e-1]);
if w<=0 & dp
   w=abs(w);
   dp=max(find([1 mod(mod(round(x*10^(dp-3*k)),10^dp),10.^(dp:-1:1))]))-1;
end
if(k)
   s=sprintf(sprintf('%%%d.%df %c',max(w-2,0),dp,f(k+7)),x*1e-3^k);
else
   s=sprintf(sprintf('%%%d.%df ',max(w-1,0),dp),x*1e-3^k);
end
