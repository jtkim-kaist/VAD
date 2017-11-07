function x=pcma2lin(p,m,s)
%PCMU2LIN Convert A-law PCM to linear X=(P,M,S)
%	lin = pcma2lin(pcma,m,s) where pcma contains a vector or matrix
%	of A-law values in the range 0 to 255.
%	No checking is performed to see that numbers are in this range.
%
%	Input values are exclusive ored with m (default=85)
%
%	Output values are divided by the scale factor s:
%
%		   s		Output Range
%
%		   1		+-4032	(integer values)
%		2017.396342	+-1.998616 (default)
%		4032		+-1
%		4096		+-0.984375 (+-1 nominal full scale)
%
%	The default value of s is 2017.396342 which equals
%	sqrt((1120^2 + 2624^2)/2). This factor follows ITU standard G.711 and
%	the sine wave with PCM-A values [225 244 244 225 97 116 116 97]
%	has a mean square value of unity corresponding to 0 dBm0.



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: pcma2lin.m 713 2011-10-16 14:45:43Z dmb $
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
  t=4.95688418E-4;
  if nargin<2 m=85; end
else
  t=1/s;
end

if m q=bitxor(p,m); else q=p; end;
k=rem(q,16);
g=floor(q/128);
e=(q-k-128*g)/16;
f=(abs(e-1)-e+1)/2;
x=(2*g-1).*(pow2(k+16.5,e)+f.*(k-15.5))*t;
