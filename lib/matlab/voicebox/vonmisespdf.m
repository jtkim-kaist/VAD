function p=vonmisespdf(x,m,k)
%VONMISESPDF Von Mises probability distribution P=(x,m,k)
%
%  Inputs:  X         matrix of input values (in radians)
%           M         mean angle of distribution (in radians)
%           K         concentration parameter
%
% Outputs:  P         matrix of probability density values (same size as X)
%                     (with no output argument, the function will plot a graph)
%
% The von Mises distribution describes the pdf of an angle over the range [0,2pi). 
% For large K, the distribution approximates a Gaussian with mean M and
% variance 1/K. For small K, the distribution is uniform.

%      Copyright (C) Mike Brookes 1997-2011
%      Version: $Id: vonmisespdf.m 4971 2014-08-05 18:34:56Z dmb $
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
if nargout>0
    p=exp(k*cos(x-m))/(2*pi*besseli(0,k));
else
    if nargin<1 || isempty(x)
        x=linspace(-pi,pi,181);
    end
    if nargin<2 || isempty(m)
        m=0;
    end
    if nargin<3 || isempty(k)
        k=[0 pow2(-1:3)];
    end
    nm=length(m);
    nk=length(k);
    np=max(nm,nk);
    pp=zeros(length(x),np);
    pl=cell(np,1);
    for i=1:np
        mi=m(1+rem(i-1,nm));
        ki=k(1+rem(i-1,nk));
        pp(:,i)=vonmisespdf(x(:),mi,ki);
        pl{i}=sprintf('\\mu=%.1f, \\kappa=%.1f',mi,ki);
    end
    plot(x,pp);
    axisenlarge([-1 -1.05]);
    legend(pl,'location','northeast');
end

