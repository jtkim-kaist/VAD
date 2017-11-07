function pf=lpcra2pf(ra,np)
%LPCAR2PF Convert AR coefs to power spectrum PF=(RA,NP)
% The routine is faster if NP+1 is a power of 2
% For RA(:,p+1) the default value of np is p and the output is PF(:,p+2)


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcra2pf.m 713 2011-10-16 14:45:43Z dmb $
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

[nf,p1]=size(ra);
if nargin<2 np=p1-1; end
pp=2*np+2;
if pp>=2*p1
   pf=abs(rfft([ra zeros(nf,pp-2*p1+1) ra(:,p1:-1:2)].').').^(-1);
else
   pf=abs(rfft([ra(:,1:np+2) ra(:,np+1:-1:2)].').').^(-1);
end   

