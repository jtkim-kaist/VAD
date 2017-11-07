function x=randiscr(p,n,a)
%RANDISCR Generate discrete random numbers with specified probabiities [X]=(P,N,A)
%
% Usage: (1) randiscr([],10)        % generate 10 uniform random binary values
%        (2) randiscr(2:6,10)       % generate 10 random numbers in the range 1:5
%                                     with probabilities [2 3 4 5 6]/20
%        (3) randiscr([],10,'abcd') % generate a string of 10 random
%                                     characters equiprobable from 'abcd'
%
% Inputs: P  vector of probabilities (not necessarily normalized) [default = uniform]
%         N  number of random values to generate [default = 1]
%         A  output alphabet [default = 1:length(p) or 0:1 if p is empty]
%
% Outputs: X  vector of not necessarily distinct values taken from alphabet A
%
% The vector P is internally normalized by dividing by its sum.
% If P is an M-dimensional matrix (and A is unspecified), then X will
% have dimensions (N,M) with the corresponding indices for each dimension.

% Somewhat similar in function to RANDSRC in the comms toolbox

%   Copyright (c) 2005-2012 Mike Brookes,  mike.brookes@ic.ac.uk
%      Version: $Id: randiscr.m 2189 2012-07-20 13:47:00Z dmb $
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
gota=nargin>2;
if nargin<1 || ~numel(p)
    if gota
        p=ones(1,length(a));
    else
        p=ones(1,2);
        a=(0:1)';
        gota=1;
    end
end
if nargin<2 || ~numel(n)
    n=1;
end
d=length(p(:)); % size of output alphabet
z=zeros(d+n-1,1); % array to hold random numbers
z(1:d)=cumsum(p(:)/sum(p(:))); % last value is actually overwritten in the next line
z(d:d+n-1)=rand(n,1);
[y,iy]=sort(z);
y(iy)=(1:d+n-1)';
m=zeros(d+n-1,1);
m(y(1:d-1))=1;
m(1)=m(1)+1;
mc=cumsum(m);
x=mc(y(d:d+n-1));
if gota
    x=a(x);
elseif length(p(:))>length(p) % need multiple dimensions
    v=x-1;
    s=cumprod(size(p));
    m=length(s);
    s(2:end)=s(1:end-1);
    s(1)=1;
    x=zeros(n,m);
    for i=m:-1:1
        x(:,i)=1+floor(v/s(i));
        v=rem(v,s(i));
    end
end