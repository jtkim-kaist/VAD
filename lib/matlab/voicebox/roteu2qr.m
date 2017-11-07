function q=roteu2qr(m,t)
%ROTEU2QR converts a sequence of Euler angles to a real unit quaternion
% Inputs:
%
%     M(1,n)   a string of n characters from the set {'x','y','z'}
%              or, equivalently, a vector whose elements are 1, 2, or 3
%     T(n,1)   n rotation angles. A positive rotation is clockwise if
%              looking along the axis away from the origin.
%
% Outputs:
%
%     Q(4,1)   output quaternion. Q is normalized to have magnitude 1 with
%              its first non-zero coefficient positive.
%
% The string M specifies the axes about which the rotations are performed.
% You cannot have the same axis in adjacent positions and so there are 12
% possibilities. Common ones are "ZXZ" and "ZYX". A positive rotation is clockwise
% if looking along the axis away from the origin; thus a rotation of +pi/2
% around Z rotates [1 0 0]' to [0 1 0]'.
% 
% Inverse conversion: If m has length 3 with adjacent characters distinct,
%                     then rotqr2eu(m,roteu2qr(m,t))=t.
%
% Inverse rotation:   qrmult(roteu2qr(m,t),roteu2qr(fliplr(m),-fliplr(t)))=+-[ 1 0 0 0]'

%
%      Copyright (C) Mike Brookes 2007-2012
%      Version: $Id: roteu2qr.m 2171 2012-07-12 07:33:03Z dmb $
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
y=[2 4 1 3 1 3 2 4; 3 2 1 4 1 4 3 2; 3 4 2 1 1 2 4 3];
% m consists of a sequence of axes e.g. 'zxy'
% and t gives the rotation angles in radians
q=[1 0 0 0]';
if ischar(m)
    m=lower(m)-'w';
end
if any(abs(m-2)>1), error('Euler axis must be x,y or z'); end
for i=1:length(m)
    x=y(m(i),:);
    b=0.5*t(i);
    c=cos(b);
    s=sin(b);
    r=zeros(4,1);
    r(x(1:2))=q(x(3:4));
    r(x(5:6))=-q(x(7:8));
    q=c*q+s*r;
end
q=q*(2*(q(find(q~=0,1))>0)-1); % force leading coefficient to be positive