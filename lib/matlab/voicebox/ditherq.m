function [y,zf]=ditherq(x,m,zi)
%DITHERQ  add dither and quantize [Y,ZF]=(X,M,ZI)
%  Inputs:
%      x   is the input signal
%	   m   specifies the mode:
%          'w'  white dither (default)
%          'h'  high-pass dither (filtered by 1 - z^-1)
%          'l'  low pass filter  (filtered by 1 + z^-1)
%          'n'  no dither

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: ditherq.m 713 2011-10-16 14:45:43Z dmb $
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

s=size(x);
n=length(x);
if nargin<3 | ~length(zi)
    zi=rand(1);
end
    if nargin<2
        m='w';
    end
if any(m=='n')
    y=round(x);
elseif any(m=='h') | any(m=='l')
    v=rand(n+1,1);
    v(1)=zi;
    zf=v(end);
    if any(m=='h')
        y=round(x(:)+v(2:end)-v(1:end-1));
    else
        y=round(x(:)+v(2:end)+v(1:end-1)-1);
    end
else
    y=round(x(:)+rand(n,2)*[1;-1]);
    zf=rand(1);                         % output a random number anyway
end
if s(1)==1
    y=y.';
end