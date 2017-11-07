function q=qrdivide(q1,q2)
%QRDIVIDE divdes two real quaternions q=[q1,q2]
%
% Inputs:
%
%     q1(4,1), q2(4,1)  Two real quaternions in the form [r, i, j, k]' where i^2=j^2=k^2=ijk=-1
%
% Outputs: 
%
%     q(4,1)   Quotient of q1/q2 such that q1=q*q2.
%              Note that q*q2 ~= q2*q since quaternion multiplication does not commute.

%      Copyright (C) Mike Brookes 2000-2008
%      Version: $Id: qrdivide.m 713 2011-10-16 14:45:43Z dmb $
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
persistent a b c d
if isempty(a)
    a=[5 8 9 10 15 13];
    b=[6 7 11 12 14 16];
    c=[1 2 3 4 6 7 11 12 16 14];
    d=[1 2 3 4 5 8 9 10 13 15];
end
if nargin<2
    %    just take the inverse of the only input argument
    q=q1/(q1'*q1);
    q(2:4)=-q(2:4);
else
    %    invert q2 and do a multiply
    q=q2/(q2'*q2);
    q(2:4)=-q(2:4);
    t=q1*q.';
    s=zeros(4,4);
    s(a)=-t(b);
    s(c)=t(d);
    q=sum(s,2);
end

