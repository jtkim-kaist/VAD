function [s,c,r,t]=atan2sc(y,x)
%ATAN2SC    sin and cosine of atan(y/x) [S,C,R,T]=(Y,X)
%
% Outputs:
%    s    sin(t) where tan(t) = y/x
%    C    cos(t) where tan(t) = y/x
%    r    sqrt(x^2 + y^2)
%    t    arctan of y/x

%      Copyright (C) Mike Brookes 2007
%      Version: $Id: atan2sc.m 713 2011-10-16 14:45:43Z dmb $
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


t=NaN;
if y == 0
    t=(x<0);
    c=1-2*t;
    s=0;
    r=abs(x);
    t=t*pi;
elseif x == 0
    s=2*(y>=0)-1;
    c=0;
    r=abs(y);
    t=s*0.5*pi;
elseif (abs(y) > abs(x))
    q = x/y;
    u = sqrt(1+q^2)*(2*(y>=0)-1);
    s = 1/u;
    c = s*q;
    r = y*u;
else
    q = y/x;
    u = sqrt(1+q^2)*(2*(x>=0)-1);
    c = 1/u;
    s = c*q;
    r = x*u;
end
if nargout>3 && isnan(t)
    t=atan2(s,c);
end