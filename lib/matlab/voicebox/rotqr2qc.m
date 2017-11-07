function qc=rotqr2qc(qr)
%ROTQR2QC converts a matrix of real quaternion vectors into complex form
%
% Inputs: 
%
%     QR(4m,...)   array of real-valued quaternions
%
% Outputs: 
%
%     QC(2m,...)   array of complex-valued quaternions
%
% The real-valued quaternion [r a b c]' becomes [r+j*b  a+j*c]'

% 
%      Copyright (C) Mike Brookes 2000-2012
%      Version: $Id: rotqr2qc.m 1616 2012-03-15 09:13:31Z dmb $
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
persistent a b
if isempty(a)
    a=[1 3 2 4];
    b=[1 1i];
end
s=size(qr);
qq=reshape(qr,4,[]);
s(1)=s(1)/2;
qc=reshape(b*reshape(qq(a,:),2,[]),s);