function x=pcmu2lin(p,s)
%PCMU2LIN Convert Mu-law PCM to linear X=(P,S)
%	lin = pcmu2lin(pcmu) where pcmu contains a vector
%	of mu-law values in the range 0 to 255.
%	No checking is performed to see that numbers are in this range.
%
%	Output values are divided by the scale factor s:
%
%		   s		Output Range
%
%		   1		+-8031	(integer values)
%		4004.2	+-2.005649 (default)
%		8031		+-1
%		8159		+-0.9843118 (+-1 nominal full scale)
%
%	The default scaling factor 4004.189931 is equal to
%	sqrt((2207^2 + 5215^2)/2) this follows ITU standard G.711.
%	The sine wave with PCM-Mu values [158 139 139 158 30 11 11 30]
%	has a mean square value of unity corresponding to 0 dBm0.



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: pcmu2lin.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<2
  t=9.98953613E-4;
else
  t=4/s;
end

m=15-rem(p,16);
q=floor(p/128);
e=(127-p-m+128*q)/16;
x=(q-0.5).*(pow2(m+16.5,e)-16.5)*t;
