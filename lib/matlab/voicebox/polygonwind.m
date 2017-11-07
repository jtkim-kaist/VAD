function w=polygonwind(p,x)
%POLYGONWIND Test if points are inside a polygon
% Inputs:
%    P(n,2)  polygon vertices
%    X(m,2)  points to test
%
% Outputs:
%    W(m)    winding number for each point
%
% For a normal polygon whose vertices are listed
% anti-clockwise, the winding number is 0 or 1 according to whether the
% point is outside or inside the polygon. The winding number will
% be -1 if the polygon vertices go colckwise around the point.


%      Copyright (C) Mike Brookes 2009
%      Version: $Id: polygonwind.m 713 2011-10-16 14:45:43Z dmb $
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
n=size(p,1);
m=size(x,1);
q=zeros(2,n+1);
q(:,1:n)=p';
q(:,n+1)=q(:,1);      % append an extra point
i=1:n;
j=2:n+1;
ym=repmat(2,m,1);
yn=repmat(2,1,n);
w=sum((2*((repmat(q(1,i).*q(2,j)-q(2,i).*q(1,j),m,1)+x(:,1)*(q(2,i)-q(2,j))+x(:,2)*(q(1,j)-q(1,i)))>0)-1).*abs((q(ym,j)>x(:,yn))-(q(ym,i)>x(:,yn))),2)/2;
if ~nargout
    w0=w==0;
    wp=w>0;
    wn=w<0;
    plot(q(1,:),q(2,:),'k-',x(w0,1),x(w0,2),'go',x(wp,1),x(wp,2),'r+',x(wn,1),x(wn,2),'bx');
    mnx=[1.05 -0.05;-0.05 1.05]*[min([p; x]);max([p; x])];
    set(gca,'xlim',mnx(:,1)','ylim',mnx(:,2)');
    title('Winding numbers: o=0, +=pos, x=neg');
end
