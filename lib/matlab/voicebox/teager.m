function y=teager(x,d,m)
%TEAGER calculate teager energy waveform Y=(X,D,M)
%
%  Inputs:  x         speech signal
%           d         dimension to apply filter along [default 1st non-singleton]
%           m         Normally Y has the same length as X and the first
%                     and last output samples are extrapolated. Setting m='x'
%                     supresses this extrapolation and Y will be two
%                     samples shorter than X
%
% Outputs:  Y         output signal: y(n)=abs(x(n))^2 - x(n+1)*conj(x(n-1))
%
% Calculates the Teager energy waveform [1]. The following waveforms give
% a constant output (independent of n) where A, B, C are real constants:
%  (a) x(n) = A*sin(B*n+C)    -->   y(n) = (A*sin(B))^2
%  (b) x(n) = A*n + B         -->   y(n) = A^2
%  (c) x(n) = A*exp(j(B*n+C)) -->   y(n) = A^2*(1-exp(2jB))
%  (d) x(n) = A*exp(B*n+C)    -->   y(n) = 0
%
% Reference:
%  [1]	J. Kaiser. On a simple algorithm to calculate the ‘energy’ of a signal.
%       In Proc IEEE Intl Conf Acoustics, Speech and Signal Processing,
%       pages 381–384, vol.1, Apr. 1990. doi: 10.1109/ICASSP.1990.115702.

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: teager.m 713 2011-10-16 14:45:43Z dmb $
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
if nargin>2 && any(m=='x')
    y=z(2:k-1,:).*conj(z(2:k-1,:))-z(3:k,:).*conj(z(1:k-2,:));
    k=k-2;              % we have lost two elements
elseif k>=4
    y=zeros(k,q);
    y(2:k-1,:)=z(2:k-1,:).*conj(z(2:k-1,:))-z(3:k,:).*conj(z(1:k-2,:));
    y(1,:)=2*y(2,:)-y(3,:);             % linearly interpolate the end points
    y(k,:)=2*y(k-1,:)-y(k-2,:);
elseif k==3
    y=repmat(x(2,:).*conj(x(2,:))-x(3,:).*conj(x(1,:)),3,1);
else
    y=zeros(k,q);
end
if d==1
    e(d)=k;
    y=reshape(y,e);
else
    r(1)=k;
    y=shiftdim(reshape(y,r),length(e)+1-d);
end