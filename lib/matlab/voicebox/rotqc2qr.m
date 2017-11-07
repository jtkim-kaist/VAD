function qr=rotqc2qr(qc)
%ROTQC2QR converts a matrix of complex quaternion row vectors into real form
%
% Inputs: 
%
%     QC(2m,n)   mxn matrix of complex-valued quaternions
%
% Outputs: 
%
%     QR(4m,n)   mxn matrix of real-valued quaternions
%
% The complex-valued quaternion [r+j*b  a+j*c] becomes [r a b c]

% 
%      Copyright (C) Mike Brookes 2000-2006
%      Version: $Id: rotqc2qr.m 713 2011-10-16 14:45:43Z dmb $
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
[m,n]=size(qc);
i=(1:2:2*m)-mod(0:m-1,2);
qr=zeros(2*m,n);
qr(i,:)=real(qc);
qr(i+2,:)=imag(qc);