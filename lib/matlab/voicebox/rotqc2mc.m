function mc=rotqc2mc(qc)
%ROTQC2MC converts a matrix of complex quaternion vectors to quaternion matrices
% Inputs: 
%
%     QC(2m,n,...)   array of complex quaternion vectors (each 2x1)
%
% Outputs: 
%
%     MC(2m,2n,...)   array of complex quaternion matrices (each 2x2)
%
% Each quarternion, r+ai+bj+ck is a 2x1 complex vector in the input array of the
% form [ r+bi ; a+ci ]. In the output array, this becomes a 2x2 complex matrix
% of the form [ r+bi -a+ci ; a+ci r-bi ].
% 
% In matrix form, quaternions can be multiplied and added using normal matrix 
% arithmetic. The complex matrix form has 4 complex elements while the real
% matrix form has 16 real elements.

% 
%      Copyright (C) Mike Brookes 2000-2012
%      Version: $Id: rotqc2mc.m 1616 2012-03-15 09:13:31Z dmb $
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
s=size(qc);
m=s(1);
qa=reshape(qc,m,[]);
n=2*size(qa,2);
mc=zeros(m,n);
ix=1:2:m;
jx=2:2:n;
mc(:,jx-1)=qa;
mc(ix,jx)=-conj(qa(ix+1,:));
mc(ix+1,jx)=conj(qa(ix,:));
s(2)=2*s(2);
mc=reshape(mc,s);