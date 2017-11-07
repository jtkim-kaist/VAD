function [r,p,q]=rotation(x,y,z)
%ROTATION Encode and decode rotation matrices
% (1) r=rotation(x,y,angle) creates a matrix that rotates vectors in the
%     plane containing x and y. A small positive angle moves x towards y.
% (2) [x,y,ang]=rotation(r) is the inverse of (1) above. The input argument r
%     is assumed to be a rotation matrix: little error checking is done.
% (3) r=rotation(axis,angle)=rotation(axis*ang) only works for a 3-dimensional
%     vector axis and creates a rotation of angle radians around axis.
% (4) [axis,ang]=rotation(r) is the inverse of (3) above for a 3x3 input matrix r

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: rotation.m 713 2011-10-16 14:45:43Z dmb $
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

l=length(x(:));
if nargin>2
   x=x(:); x=x/sqrt(x'*x);
   y=y(:); y=y-y'*x*x; y=y/sqrt(y'*y);
   r=eye(l)+(cos(z)-1)*(x*x'+y*y')+sin(z)*(y*x'-x*y');
elseif l==3
   x=x(:);
   lx=sqrt(x'*x);
   if nargin==1
      y=lx;
   end
   x=x/lx;
   xx=x*x';
   s=zeros(3,3);
   s([6 7 2])=x;
   s([8 3 4])=-x;
   r=xx+cos(y)*(eye(3)-xx)+sin(y)*s;
else
   [v,d]=eig(x);
   [e,j]=sort(real(diag(d)));
   j1=j(1);
   an=angle(d(j1,j1));
   q=an;
   sq=sqrt(2);
   r=imag(v(:,j1))*sq;
   if r==0
      p=v(:,j1);
      r=v(:,j(2));
   else
      p=real(v(:,j1))*sq;
   end
   if nargout==2
      r=cross(r,p);
      p=an;
   end
end