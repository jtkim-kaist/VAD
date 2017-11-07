function p=normcdflog(x,m,s)
%NORMCDFLOG calculates log of Normal Cumulative Distribution function p=(x,m,s)
%
% Inputs: 
%
%     X        Input data (vector or matrix)
%     M        Mean of Normal distribution [default 0]
%     S        Std deviation of Normal distribution [default 1]
%
% Outputs: 
%
%     P        P = log(normcdf(X)); same size as X
%
% The routine gives accurate values even if X is large and negative

%      Copyright (C) Mike Brookes 2016
%      Version: $Id: normcdflog.m 8839 2016-10-19 15:03:19Z dmb $
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
persistent a b
if isempty(a)
    a=0.996;
%   f=@(x) -0.5*(x.^2+log(2*pi))-real(log(-(a+x.^2)./x))-real(log(normcdf(x)));
%   b=fzero(f,-22);         % cutoff value for conventional formula
    b=-22.2491306156561;    % precalculated value
end
if nargin>2
    x=(x-m)/s;
elseif nargin>1
    x=x-m;
end
t=x<b;                                                  % mask for large negative values
p=zeros(size(x));
p(~t)=real(log(0.5*erfc(-x(~t).*sqrt(0.5))));           % use this formula normally
p(t)=-0.5*(x(t).^2+log(2*pi))-real(log(-x(t)-a./x(t))); % use this approximation for large negative x

