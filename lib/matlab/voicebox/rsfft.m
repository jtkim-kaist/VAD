function x=rsfft(y,n)
%RSFFT    fft of a real symmetric spectrum X=(Y,N)
% Y is the "first half" of a symmetric real input signal and X is the
% "first half" of the symmetric real fourier transform.
% If the length, N, of the full signal is even, then the "first half"
% contains 1+N/2 elements (the first and last are excluded from the reflection).
% If N is odd, the "first half" conatins 0.5+N/2 elements and only the first
% is excluded from the reflection.
% If N is specified explicitly, then Y will be truncated of zero-padded accordingly.
% If N is omitted it will be taken to be 2*(length(Y)-1) and is always even.
%
% If Y is a matrix, the transform is performed along each column
%
% The inverse function is y=rsfft(x,n)/n

% Could be made faster for even n by using symmetry

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: rsfft.m 5024 2014-08-22 17:05:55Z dmb $
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

if ~isreal(y) error('RSFFT: Input must be real'); end
fl=size(y,1)==1;
if fl y=y(:); end
[m,k]=size(y);
if nargin<2 n=2*m-2;
else
    mm=1+fix(n/2);
    if mm>m y=[y; zeros(mm-m,k)];
    elseif mm<m y(mm+1:m,:)=[];
    end
    m=mm;
end
x=real(fft([y;y(n-m+1:-1:2,:)]));
x(m+1:end,:)=[];
if fl x=x.'; end
