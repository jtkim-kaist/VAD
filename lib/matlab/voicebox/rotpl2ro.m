function r=rotpl2ro(u,v,t)
%ROTPL2RO find matrix to rotate in the plane containing u and v r=[u,v,t]
% Inputs:
%
%     U(n,1) and V(n,1) define a plane in n-dimensional space
%     T is the rotation angle in radians from U towards V. If T
%       is omitted it will default to the angle between U and V
%
% Outputs:
%
%     R(n,n)   Rotation matrix

%
%      Copyright (C) Mike Brookes 2007
%      Version: $Id: rotpl2ro.m 713 2011-10-16 14:45:43Z dmb $
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

u=u(:);
    n=length(u);
v=v(:);
l=sqrt(u'*u);
if l==0, error('input u is a zero vector'); end
u=u/l;      % normalize
q=v-v'*u*u;        % q is orthogonal to x
l=sqrt(q'*q);
if l==0          % u and v are colinear or v=zero
    [m,i]=max(abs(u));
    q=zeros(n,1);
    q(1+mod(i(1),n))=1;  % choose next available dimension
    q=q-q'*u*u;  % q is orthogonal to x
    l=sqrt(q'*q);
end
q=q/l;          % normalize
if nargin<3
    [s,c]=atan2sc(v'*q,v'*u);
    r=eye(n)+(c-1)*(u*u'+q*q')+s*(q*u'-u*q');
else
    r=eye(n)+(cos(t)-1)*(u*u'+q*q')+sin(t)*(q*u'-u*q');
end

