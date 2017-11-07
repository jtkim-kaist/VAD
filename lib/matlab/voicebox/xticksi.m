function s=xticksi(ah)
%XTIXKSI labels the x-axis of a plot using SI multipliers S=(AH)
%
%  Inputs:  AH       axis handle [default: current axes]
%
% Outputs:  S        optional global SI multiplier (see usage below)
%
% Usage:   (1) plot(...);
%              xticksi;
%
%          (2) plot(...);
%              xlabel(['Frequency (' xticksi 'Hz)']);
%
% The first form will label the tick marks on the x-axis of the current plot
% using SI multipliers where appropriate. This is particularly useful for log
% plots which MATLAB does not label very well.
% The second form will, if possible, use a single SI multiplier for all the tick
% marks; this global multiplier can be incorporated into the axis label as shown.

%	   Copyright (C) Mike Brookes 2009
%      Version: $Id: xticksi.m 713 2011-10-16 14:45:43Z dmb $
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
if ~nargin
    ah=gca;
end
if nargout
s=xyzticksi(1,ah);
else
    xyzticksi(1,ah);
end