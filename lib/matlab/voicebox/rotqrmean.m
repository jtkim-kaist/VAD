function [y,s,v]=rotqrmean(q)
%ROTQRMEAN calculates the mean rotation of a quaternion array [y,s]=[q]
%
% Inputs:   q(4,n)    normalized real quaternion array
%
% Outputs:  y(4,1)    normalized mean quaternion
%           s(1,n)    sign vector such that y=q*s', y=y/sqrt(y.'*y)
%           v         average squared deviation from the mean quaternion
%
% Since quaternions represent a rotation only to within a sign ambiguity
% we need to select +1 or -1 for each one when calculating the mean.
% This routine selects the sign for each quaternion to maximize the norm
% of their sum or, equivalently, to minimize their variance.
%
%      Copyright (C) Mike Brookes 2011-2012
%      Version: $Id: rotqrmean.m 2188 2012-07-20 13:46:29Z dmb $
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
mmax=10;                % number of n-best hypotheses to keep
nq=size(q,2);
mkx=zeros(nq,mmax);     % save signs: 0=+, 1 = -
mprev=ones(nq,mmax);    % save back pointers
msum=zeros(4,2*mmax);
msum(:,1)=q(:,1);       % current values of sum
ix=1:mmax;
jx=mmax+1:2*mmax;
r=ones(1,mmax);
for i=2:nq
    msum(:,jx)=msum(:,ix)-q(:,i(r));
    msum(:,ix)=msum(:,ix)+q(:,i(r));
    [vx,kx]=sort(sum(msum.^2,1),2,'descend');
    mkx(i,:)=kx(ix);    % negative is > mmax
    msum(:,ix)=msum(:,kx(ix));  % save mmax sums having highest norms
end
y=msum(:,1);            % unnormalized mean
y=y/sqrt(y.'*y);
if nargout>1            % do traceback
    s=zeros(1,nq);
    k=1;
    for i=nq:-1:2
        s(i)=(mkx(i,k)>mmax);
        k=mkx(i,k)-mmax*s(i);
    end
    s=1-2*s;
    v=sum(mean((q-y*s).^2,2));
end
   
    