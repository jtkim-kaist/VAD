function y=skew3d(x,m)
%SKEW3D Convert between a vector and the corresponding skew-symmetric matrix
%
% Inputs:   x   input vector or matrix
%                size(x) must equal [3 1], [3 3], [6 1] or [4 4]
%           m   m string:
%               'n'  normalize the vector to have unit magnitude
%               'z'  orthoganlize the vector so that x'Jy=0
%
% Outputs:  y   output matrix or vector
%                size(y) = [3 3], [3 1], [4 4] or [6 1] respectively
%                Note that skew3d() is its own inverse: skew3d(skew3d(x)) =  x
%
% 3D Euclidean space
% ------------------
%    If A and B are 3x1 vectors then the vector cross product is given by
%    skew3d(A)*B = cross(A,B) = A x B. This relationship is widely used
%    in computer vision.
%
% 3D Projective space
% -------------------
% In 3D projective space, a line has 4 degrees of freedom and may be
% represented by its homogeneous 6x1 Plucker vector, A, or 4x4 Plucker
% matrix B=skew3d(A).
% The 6x1 Plucker vector loses one degree of freedom because it is
% homogeneous (i.e. independent of a non-zero scale factor) and another
% because it must satisfy A'*flipud(A)=0. Setting the 'n' and 'z' options
% in the second input parameter will remove these redundancies by forcing
% A'*A=1 and A'*flipud(A)=0.

%      Copyright (C) Mike Brookes 1998-2012
%      Version: $Id: skew3d.m 2152 2012-07-06 15:41:30Z dmb $
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

[j,k]=size(x);
mn=nargin>1 && any(m=='n');  % normalize
mz=nargin>1 && any(m=='z');  % orthoganalize
if j==3
    if k==1
        if mn && x'*x>0
            x=x/sqrt(x'*x);
        end
        y=zeros(3,3);
        y([6 7 2])=x(:)';
        y([8 3 4])=-x(:)';
    elseif k==3
        y=x([6 7 2]');
        if mn && y'*y>0
            y=y/sqrt(y'*y);
        end
    else
        error('size(x) must be [3 1], [3 3], [6 1] or [4 4]');
    end
elseif j==6 && k==1
    x=x(:);
    u=x(1:3);
    v=x(6:-1:4);
    if mz && u'*u>0 && v'*v>0  % orthoganalize
        v = v - (u'*v)/(2*u'*u)*u;
        x = [u-(v'*u)/(v'*v)*v; v([3 2 1])];
    end
    if mn && x'*x>0
        x=x/sqrt(x'*x);
    end
    y=zeros(4,4);
    y([5 9 13 10 8 15])=x(:)';
    y([2 3 4 7 14 12])=-x(:)';
elseif j==4 && k==4
    u=x([5 9 13]');
    v=x([15 8 10]');
    if mz && u'*u>0 && v'*v>0  % orthoganalize
        v = v - (u'*v)/(2*u'*u)*u;
        y = [u-(v'*u)/(v'*v)*v; v([3 2 1])];
    else
        y = [u; v([3 2 1])];
    end
    if mn && y'*y>0
        y=y/sqrt(y'*y);
    end
else
    error('size(x) must be [3 1], [3 3], [6 1] or [4 4]');
end
