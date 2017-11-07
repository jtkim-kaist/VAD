function [y,t]=schmitt(x,thresh,minwid)
% Pass input signal X through a schmitt trigger
% SCHMITT(X,[LOW HIGH]) gives low and high thresholds. LOW and HIGH can be
% scalars or can be vectors specifiying different thresholds for each X element.
%
% SCHMITT(X,HYSTERESIS) specifies the thresholds as MAX-DELTA and MIN+DELTA where
% DELTA=(MAX-MIN)*(1-HYSTERESIS)/2 and MAX and MIN are the max and min of X.
% HYSTERESIS must be in the range 0 to 1 and represents the fraction of MAX-MIN between
% the two threshold values; default is 0.5.
%
% MINWID specifies the minimum width of a pulse (in samples). Pulses thinner than this will
% be ignored.
%
% Output Y takes values -1, +1 according to whether X<LOW or X>HIGH most recently.
% Y may be 0 for an initial segment if neither condition is initially true.
% For [Y,T]=SCHMITT(...) Y contains alternate +1 and -1 values and T contains
% the sample numbers at which X crossed the thresholds.


%      Copyright (C) Mike Brookes 1998
%      Version: $Id: schmitt.m 713 2011-10-16 14:45:43Z dmb $
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
    minwid=0;
    if nargin<2
        thresh=0.5;
    end
end
if length(thresh)<2
    xmax=max(x);
    xmin=min(x);
    low=(xmax-xmin)*(1-thresh)/2;
    high=xmax-low;
    low=xmin+low;
else
    low=thresh(1);
    high=thresh(2);
end
c=(x>high)-(x<low);
c(2:end)=c(2:end).*(c(2:end)~=c(1:end-1));
t=find(c);
t(1+find(c(t(2:end))==c(t(1:end-1))))=[]; % remove duplicates
if minwid>=1
    t(t(2:end)-t(1:end-1)<minwid)=[];
    t(1+find(c(t(2:end))==c(t(1:end-1))))=[]; % remove duplicates
end
if nargout>1 y=c(t);
else
    y=zeros(size(c));
    if ~isempty(t)
        y(t)=2*c(t);
        y(t(1))=c(t(1));
        y=cumsum(y);
    end
end
if ~nargout
        xmax=max(x);
    xmin=min(x);
    if high-low<0.1*(xmax-xmin)
        high=xmax;
        low=xmin
    end
    plot([x(:) low+(y(:)+1)*(high-low)/2]);
end