function q=qrdotmult(q1,q2)
%QRDOTMULT multiplies together two real quaternions arrays q=[q1,q2]
%
% Inputs:  q1(4n,...)
%          q2(4n,...)   Two real quaternion arrays of equal size
%
% Outputs: q(4n,...)    The Hadamard product (i.e. .*) of the input arrays
%
%      Copyright (C) Mike Brookes 2000-2012
%      Version: $Id: qrdotmult.m 1618 2012-03-15 09:14:33Z dmb $
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
persistent a b c
if isempty(a)
    a=[1 2 3 4 2 1 4 3 3 4 1 2 4 3 2 1];
    b=[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
    c=[2 3 4 7 12 14];
end
s=size(q1);
qa=reshape(q1,4,[]);
qb=reshape(q2,4,[]);
q=qa(a,:).*qb(b,:);
q(c,:)=-q(c,:);
q=reshape(sum(reshape(q,4,[]),1),s);
