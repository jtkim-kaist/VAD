function q=qrdotdiv(x,y)
%QRDOTDIV divides two real quaternions arrays elementwise q=[x,y]
%
% Inputs:  x(4n,...)
%          y(4n,...)   Two real quaternion arrays of equal size
%
% Outputs: q(4n,...)    The Hadamard quaotient (i.e. ./) of the input arrays
%                       If y is omitted then q = x^(-1)
%
%      Copyright (C) Mike Brookes 2000-2012
%      Version: $Id: qrdotdiv.m 1689 2012-03-22 21:45:41Z dmb $
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
    c=[6 8 10 11 15 16];
end
s=size(x);
x=reshape(x,4,[]);
if nargin<2 % one argument - just invert x
    m=sum(x.^2,1);
    x=x./m(ones(4,1),:);
    x(2:4,:)=-x(2:4,:);
    q=reshape(x,s);
else  % multiply by conjugate of y and then divide by |y|^2
    y=reshape(y,4,[]);
    m=sum(y.^2,1);
    q=x(a,:).*y(b,:);
    q(c,:)=-q(c,:);
    q=reshape(sum(reshape(q,4,[]),1),s)./m(ones(4,1),:);;
end
