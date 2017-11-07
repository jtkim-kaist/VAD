function m = rnsubset(k,n)
%RNSUBSET choose k distinct random integers from 1:n M=(K,N)
%
%  Inputs:
%
%    K is number of disinct integers required from the range 1:N
%    N specifies the range - we must have K<=N
%
%  Outputs:
%
%    M(1,K) contains the output numbers

%      Copyright (C) Mike Brookes 2006
%      Version: $Id: rnsubset.m 713 2011-10-16 14:45:43Z dmb $
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
if k>n
    error('rnsubset: k must be <= n');
end
% We use two algorithms according to the values of k and n
[f,e]=log2(n);
if k>0.03*n*(e-1)
[v,m]=sort(rand(1,n)); % for large k, just do a random permutation
else
    v=ceil(rand(1,k).*(n:-1:n-k+1));
    m=1:n;
    for i=1:k
        j=v(i)+i-1;
        x=m(i);
        m(i)=m(j);
        m(j)=x;
    end
end
m=m(1:k);
