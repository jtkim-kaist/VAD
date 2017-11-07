function e=rotro2eu(m,ro)
%ROTRO2EQ converts a 3x3 rotation matrix into the corresponding euler angles
% Inputs: 
%
%     M(1,3)   a string of 3 characters from the set {'x','y','z'} e.g. "zxz" or "zyx"
%              or, equivalently, a vector whose elements are 1, 2 or 3
%     RO(3,3)  3x3 rotation matrix
%
% Outputs:
%
%     E(3,1)   3 euler angles in the range +-pi. A positive rotation is clockwise
%               if looking along the axis away from the origin.
%
% The string M specifies the axes (fixed in space) about which the rotations of
% an object are performed. You cannot have the same axis in adjacent positions
% and so there are 12 possibilities. Common ones are "ZXZ" and "ZYX".
% If you want the axes to move with the object, you should reverse the
% ordering of both "m" and "e".
%
% There is some reduncancy in euler angle values:
%  (i)  If m(1)==m(3) then e=[a b c] and e=[a+-pi -b c+-pi] are equivalent.
%       The output of this routine will always have b>=0;
%  (ii) If m(1)~=m(3) then e=[a b c] and e=[a+-pi pi-b c+-pi] are equivalent.
%       The output of this routine will always have |b|<=pi/2

% 
%      Copyright (C) Mike Brookes 2007
%      Version: $Id: rotro2eu.m 2187 2012-07-20 13:45:35Z dmb $
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

e=zeros(3,1);
if ischar(m)
    m=lower(m)-'w';
end
if numel(m)~=3 | any(abs(m-2)>1), error('Euler axis must be an x,y or z triplet'); end
u=m(1);
v=m(2);
w=m(3);
if sum(m==v)>1, error('Consecutive Euler axes must differ'); end
% first we rotate around w to null element (v,u) with respect to element (!vw,u) of rotation matrix
g=2*mod(u-v,3)-3;   % +1 if v follows u or -1 if u follows v
h=2*mod(v-w,3)-3;   % +1 if w follows v or -1 if v follows w
[s,c,r,e(3)]=atan2sc(h*ro(v,u),ro(6-v-w,u));
r2=ro;
ix=1+mod(w+(0:1),3);
r2(ix,:)=[c s; -s c]*ro(ix,:);
% next we rotate around v to null element (!uv,u) with repect to element (u,u) of rotation matrix
e(2)=atan2(-g*r2(6-u-v,u),r2(u,u));
% finally we rotate around u to null element (v,!uv) with respect to element (!uv,!uv) = element (v,v)
e(1)=atan2(-g*r2(v,6-u-v),r2(v,v));
if (u==w && e(2)<0) || (u~=w && abs(e(2))>pi/2)  % remove redundancy
    mk=u~=w;
    e(2)=(2*mk-1)*e(2);
    e=e-((2*(e>0)-1) .* [1; mk; 1])*pi;
end

