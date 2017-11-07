function [t,s]=zerocros(y,m,x)
%ZEROCROS finds the zeros crossings in a signal [T,S]=(Y,M,X)
% Inputs:  y = input waveform
%          m = mode string containing:
%              'p' - positive crossings only
%              'n' - negative crossings only
%              'b' - both (default)
%              'r' - round to sample values
%          x = x-axis values corresponding to y [default 1:length(y)]
%
% Outputs: t = x-axis positions of zero crossings
%          s = estimated slope of y at the zero crossing
%
% This routine uses linear interpolation to estimate the position of a zero crossing
% A zero crossing occurs between y(n) and y(n+1) iff (y(n)>=0) ~= (y(n+1)>=0)

% Example: y=sin(2*pi*(0:1000)/200); y(1:100:1001)=0; zerocros(y);
% Note that we get a zero crossing at the end but not at the start.

%	   Copyright (C) Mike Brookes 2003-2015
%      Version: $Id: zerocros.m 6563 2015-08-16 16:56:24Z dmb $
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

if nargin<2 || ~numel(m)
    m='b';
end
s=y>=0;
k=s(2:end)-s(1:end-1);
if any(m=='p')
    f=find(k>0);
elseif any(m=='n')
    f=find(k<0);
else
    f=find(k~=0);
end
s=y(f+1)-y(f);
t=f-y(f)./s;
if any(m=='r')
    t=round(t);
end
if nargin>2
    tf=t-f; % fractional sample
    t=x(f).*(1-tf)+x(f+1).*tf;
    s=s./(x(f+1)-x(f));
end
if ~nargout
    n=length(y);
    plot(1:n,y,'-',t,zeros(length(t),1),'o');
end
