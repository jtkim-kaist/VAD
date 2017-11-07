function p=lin2pcma(x,m,s)
%LIN2PCMA Convert linear PCM to A-law P=(X,M,S)
%	pcma = lin2pcma(lin) where lin contains a vector
%	or matrix of signal values.
%	The input values will be converted to integer
%	A-law pcm vlues in the range 0 to 255 and the XORed with m
%	(default m=85).
%	
%	Input values are multiplied by the scale factor s:
%
%		   s		Input Range
%
%		   1		+-4096
%		2017.396342	+-2.03033976 (default)
%		4096		+-1
%
%	Input values outside the selected range will be clipped.
%
%	The default value of the scale factor is 2017.396342 which equals
%	sqrt((1120^2 + 2624^2)/2). This factor follows ITU standard G.711 and
%	the sine wave with PCM-A values [225 244 244 225 97 116 116 97]
%	has a mean square value of unity corresponding to 0 dBm0.
%
%	See also PCMA2LIN, LIN2PCMA, LIN2PCMU



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: lin2pcma.m 713 2011-10-16 14:45:43Z dmb $
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
  s=2017.396342;
  if nargin<2 m=85; end
end

y=x*pow2(s,-6);
y=(abs(y+63)-abs(y-63))/2;
q=floor((y+64)/64);
[a,e]=log2(abs(y));
d = (e+abs(e))/2;
p=128*q+16*d+floor(pow2(a,e-d+5));
if m p=bitxor(p,m); end;
