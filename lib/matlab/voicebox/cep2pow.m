function [m,c]=cep2pow(u,v,mode)
%CEP2POW convert cepstral means and variances to the power domain
% Inputs:
%    u: vector giving the cepstral means with u(1) the 0'th cepstral coefficient
%    v: cepstral covariance matrix or else a vector containing the diagonal elements 
% mode: 'c'  pow=exp(irdct(cep))    [default]
%       'f'  pow=exp(rsfft(cep)/n)  [fft length even]
%       'fo' pow=exp(rsfft(cep)/n)  [fft length odd]
%       'i'  pow=exp(cep)           [ no transformation ]
%
% Outputs:
%    m: row vector giving means in the power domain
%    c: covariance matrix in the power domain

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: cep2pow.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<3 mode='c'; end
if min(size(v))==1
   v=diag(v);
end
u=u(:)';    % force u to be a row vector
if any(mode=='f')
   n=2*length(u)-2;
   if any(mode=='o')
      n=n+1;
   end
   p=rsfft(u',n)/n;
   q=rsfft(rsfft(v,n)',n)/n^2;
elseif any(mode=='i')
    p=u';
    q=v';
else
   p=irdct(u');
   q=irdct(irdct(v)');
end
m=exp(p+0.5*diag(q))';
c=(m'*m).*(exp(q)-1);