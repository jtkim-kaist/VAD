function p=lognmpdf(x,m,v)
%LOGNMPDF calculate pdf of a multivariate lognormal distribution P=(X,M,V)
%
%  Inputs:  X(N,D)   are the points at which to calculate the pdf (one point per row)
%           M(D)     is the mean vector of the distribution [default M = ones]
%           V(D,D)   is the covariance matrix of the distribution. If V is diagonal
%                    it may be given as a vector [default V = identity matrix]
%
% Outputs:  P(N,1)   is the pdf at each row of X
%
% Example: lognmpdf(linspace(0,10,1000)',2);

%	   Copyright (C) Mike Brookes 1995
%      Version: $Id: lognmpdf.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<3
    if nargin<2
        m=ones(1,size(x,2));
    end
    v=eye(length(m));
end
if(size(x,2)~=length(m)) | (size(x,2)~=length(v))
    error('Number of columns must match mean and variance dimensions');
end
[u,k]=pow2cep(m,v,'i'); % convert to log domain
p=zeros(size(x,1),1);
c=prod(x,2);
q=c>0;
p(q)=mvnpdf(log(x(q,:)),u,k)./c(q);

if ~nargout & (length(u)==1)
    plot(x,p);
end