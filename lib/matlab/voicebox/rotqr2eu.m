function e=rotqr2eu(m,q)
%ROTQR2EQ converts a real unit quaternion into the corresponding euler angles
% Inputs: 
%
%     M(1,3)   a string of 3 characters from the set {'x','y','z'}
%              or, equivalently, a vector whose elements are all 1, 2 or 3
%     Q(3,3)   3x3 rotation matrix
%
% Outputs:
%
%     E(3,1)   3 euler angles
%
% The string M specifies the axes about which the rotations are performed.
% You cannot have the same axis in adjacent positions and so there are 12
% possibilities. Common ones are "ZXZ" and "ZYX".

% Suggestions:
%   (1) Might be quicker to convert to a matrix and do it in that domain

% 
%      Copyright (C) Mike Brookes 2007
%      Version: $Id: rotqr2eu.m 713 2011-10-16 14:45:43Z dmb $
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

% test: ea=2*pi*(rand(3,1)-0.5);m='yxy';q=roteu2qr(m,ea); e=rotqr2eu(m,q);[[ea*180/pi; 0] [e*180/pi; 0] q roteu2qr(m,e)]
% redundancy: If m(1)==m(3) then q= [a b c] and [a+-pi -b c+-pi] are equivalent. The output always has b>=0
%             If m(1) ~=m(3) then q=[a b c] and [a+-pi pi-b c+-pi] are equivalent. The output always has |b|<=pi/2
% all angles are in the range +-pi
e=zeros(3,1);
y=[2 4 1 3 1 3 2 4; 3 2 1 4 1 4 3 2; 3 4 2 1 1 2 4 3];
if ischar(m)
    m=lower(m)-'w';
end
if any(abs(m-2)>1), error('Euler axis must be x,y or z'); end
u=m(1)+1;
v=m(2)+1;
w=m(3)+1;
% first we rotate around w to null element (v,u) with respect to element (!vw,u) of rotation matrix
g=2*mod(u-v,3)-3;
ss=(2*mod(v-w,3)-3)*(q(v)*q(u)+g*q(9-u-v)*q(1));
if u==w                 % if u==w then (!vw,u) is off-diagonal
    cc=q(9-v-w)*q(u)-g*q(v+w-u)*q(1);
else               % if u~=w then (!vw,u)=(u,u) is on diagonal
    cc=q(1)^2+q(u)^2-0.5;
end
[ss,cc,rr,t]=atan2sc(ss,cc);
if cc>0
    c=sqrt(0.5*(1+cc));
    s=0.5*ss/c;
else
    s=sqrt(0.5*(1-cc));
    c=0.5*ss/s;
end
% s=sin(t/2);
% c=cos(t/2);
x=y(w-1,:);
r=zeros(4,1);
r(x(1:2))=q(x(3:4));
r(x(5:6))=-q(x(7:8));
q2=c*q-s*r;
% next we rotate around v to null element (!uv,u) with repect to element (u,u) of rotation matrix
ss2=-g*q2(9-u-v)*q2(u)+q2(v)*q2(1);     % always off-diagonal
cc2=q2(1)^2+q2(u)^2-0.5;     % always on-diagonal
[s2,c2,rr,t2]=atan2sc(ss2,cc2);
x2=y(v-1,:);
r2=zeros(4,1);
r2(x2(1:2))=q2(x2(3:4));
r2(x2(5:6))=-q2(x2(7:8));
q3=c2*q2-s2*r2;
if q3(1)<0, q3=-q3; end
t3=2*atan2(q3(u),q3(1));
e(1)=t3;
e(2)=t2;
e(3)=t;
if (u==w && t2<0) || (u~=w && abs(t2)>pi/2)  % remove redundancy
    mk=u~=w;
    e(2)=(2*mk-1)*t2;
    e=e-((2*(e>0)-1) .* [1; mk; 1])*pi;
end

