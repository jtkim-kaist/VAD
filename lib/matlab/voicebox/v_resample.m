function [y,h]=v_resample(x,p,q,n,b)
%V_RESAMPLE Resample and remove end transients [y,h]=(x,p,q,n,b)
%
% This multiplies the sample rate of x by p/q. It is identical to resample()
% except that the initial and final filter transients are removed.
% The number of ouput samples will be length(x)*ceil(p/q) - 2*ceil(n*max(p/q,1))
% where the filter length n has a default value of 10.
%
% Inputs:  x    input signal (or multiple signals with one per column)
%          p,q  sampling rate is multiplied by p/q (p and q must be +ve integers)
%          n    length of filter [default: 10]
%          b    Kaiser window parameter beta [default: 5]
%
% Outputs: y    resampled output signal
%          h    filter used (at a rate of p times the input sample rate)

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: v_resample.m 5119 2014-09-11 07:22:12Z dmb $
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
if p==q
    y=x;
    h=1;
else
    if nargin < 5,  b = 5;  end   % design parameter for Kaiser window LPF
    if nargin < 4,   n = 10;   end
    [y,h]=resample(x,p,q,n,b);
    m=ceil(n*max(p/q,1));
    if size(x,1)==1 % x is a row vector
        y=y(m+1:end-m);
    else
        y=y(m+1:end-m,:);
    end
end