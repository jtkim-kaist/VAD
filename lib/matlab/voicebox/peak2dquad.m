function [v,xy,t,m]=peak2dquad(z)
%PEAK2DQUAD find quadratically-interpolated peak in a 2D array
%
% Note: This routine has been superceeded by quadpeak
%
%  Inputs:  Z(m,n)   is the input array
%
% Outputs:  V        is the peak value
%           XY(2,1)  is the position of the peak (in fractional subscript values)
%           T        is -1, 0, +1 for maximum, saddle point or minimum
%           M        defines the fitted quadratic: z = [x y 1]*M*[x y 1]'

%	   Copyright (C) Mike Brookes 2008
%      Version: $Id: peak2dquad.m 713 2011-10-16 14:45:43Z dmb $
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

[v,xy,t,m]=quadpeak(z);