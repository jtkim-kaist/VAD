function u=glotros(d,t,p)
%GLOTROS  Rosenberg glottal model U=(D,T,P)
% d is derivative of flow waveform
% t is in fractions of a cycle
% p has parameters
%	p(1)=closure time
%	p(2)=+ve/-ve slope ratio



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: glotros.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin < 2
  tt=(0:99)'/100;
else
  tt=mod(t,1);
end
u=zeros(size(tt));
de=[0.6 0.5]';
if nargin < 3
  p=de;
elseif length(p)<2
  p=[p(:); de(length(p)+1:2)];
end
pp=p(1)/(1+p(2));
ta=tt<pp;
tb=tt<p(1) & ~ta;
wa=pi/pp;
wb=0.5*pi/(p(1)-pp);
fb=wb*pp;
if d==0
  u(ta)=0.5*(1-cos(wa*tt(ta)));
  u(tb)=cos(wb*tt(tb)-fb);
elseif d==1
  u(ta)=0.5*wa*sin(wa*tt(ta));
  u(tb)=-wb*sin(wb*tt(tb)-fb);
elseif d==2
  u(ta)=0.5*wa^2.*cos(wa*tt(ta));
  u(tb)=-wb^2*cos(wb*tt(tb)-fb);
else
  error('Derivative must be 0,1 or 2');
end

