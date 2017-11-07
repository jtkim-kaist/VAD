function [v,y,w]=nearnonz(x,d)
%NEARNONZ replace each zero element with the nearest non-zero element [V,Y,W]=nearnonz(X,D)
%
%  Inputs:  x         input vector, matrix or larger array
%           d         dimension to apply filter along [default 1st non-singleton]
%
% Outputs:  v         v is the same size as x but with each zero entry replaced by
%                     the nearest non-zero value along dimension d
%                     elements equidistant from two non-zero entries will be taken
%                     from the higher index
%           y         y is the same size as x and gives the index along dimension d
%                     from which the corresponding entry in v was taken
%                     If there are no non-zero entries, then the corresponding
%                     elements of y will be zero.
%           w         w is the same size as x and gives the distance (+ or -) to the
%                     nearest non-zero entry in x

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: nearnonz.m 713 2011-10-16 14:45:43Z dmb $
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
e=size(x);
p=prod(e);
if nargin<2             % if no dimension given, find the first non-singleton
    d=find(e>1,1);
    if ~numel(d)
        d=1;
    end
end
k=e(d);                 % size of active dimension
q=p/k;                  % size of remainder
if d==1
    z=reshape(x,k,q);
else
    z=shiftdim(x,d-1);
    r=size(z);
    z=reshape(z,k,q);
end
xx=z~=0;
cx=cumsum(xx);
[i,j]=find(z);
qq=cx(xx);
pos=full(sparse(qq,j,i,k,q)); % list the positions of non-zero elements in each column
mp=ceil((pos(1:end-1,:)+pos(2:end,:))*0.5); % find the mid point between consecutive non-zero elements
[i2,j2]=find(pos(2:end,:)>0);
zz=1+cumsum(full(sparse(mp(pos(2:end,:)>0),j2,1,k,q)));
y=pos(zz+repmat((0:q-1)*k,k,1));
v=z(max(y,1)+repmat((0:q-1)*k,k,1));
w=y-repmat((1:k)',1,q);
w(y==0)=0;
if d==1
    y=reshape(y,e);
    v=reshape(v,e);
    w=reshape(w,e);
else
    y=shiftdim(reshape(y,r),length(e)+1-d);
    v=shiftdim(reshape(v,r),length(e)+1-d);
    w=shiftdim(reshape(w,r),length(e)+1-d);
end
