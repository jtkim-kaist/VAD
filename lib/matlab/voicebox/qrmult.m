function q=qrmult(q1,q2)
%QRMULT multiplies together two real quaternions matrices q=[q1,q2]
%
% Inputs:   q1(4m,n)  Two real quaternions arrays. Either array can
%           q2(4n,r)  also be a scalar quaternion.
%
% Outputs:   q(4m,r)  Matrix product of q1 and q2

%      Copyright (C) Mike Brookes 2000-2012
%      Version: $Id: qrmult.m 1617 2012-03-15 09:14:01Z dmb $
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
s1=size(q1);
s2=size(q2);
if isequal(s1,[4 1])
    q=qrdotmult(repmat(q1,s2(1)/4,s2(2)),q2);
elseif isequal(s2,[4 1])
    q=qrdotmult(q1,repmat(q2,s1(1)/4,s1(2)));
else
    q=rotqr2mr(q1)*q2;
end
