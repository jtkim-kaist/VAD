function qr=rotmr2qr(mr)
%ROTMR2QR converts a matrix of real quaternion matrices to quaternion vectors
% Inputs: 
%
%     MR(4m,4n)   mxn matrix of real quaternion matrices (each 4x4)
%
% Outputs: 
%
%     QR(4m,n)   mxn matrix of real quaternion vectors (each 4x1)
%
% In matrix form, quaternions can be multiplied and added using normal matrix 
% arithmetic. Each element of an mxn matrix of quaternions is itself a 4x4 block
% so the total dimension of MR is 4m x 4n.

% 
%      Copyright (C) Mike Brookes 2000-2006
%      Version: $Id: rotmr2qr.m 713 2011-10-16 14:45:43Z dmb $
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
qr=mr(:,1:4:end);