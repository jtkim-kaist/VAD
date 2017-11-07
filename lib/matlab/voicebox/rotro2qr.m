function q=rotro2qr(r)
%ROTRO2QR converts a 3x3 rotation matrix to a real quaternion
% Inputs: 
%
%     R(3,3)   Input rotation matrix
%
% Outputs: 
%
%     Q(4,1)   normalized real-valued quaternion
%
% In the quaternion representation of a rotation, and q(1) = cos(t/2) 
% where t is the angle of rotation in the range 0 to 2pi
% and q(2:4)/sin(t/2) is a unit vector lying along the axis of rotation
% a positive rotation about [0 0 1] takes the X axis towards the Y axis.

%      Copyright (C) Mike Brookes 2007
%      Version: $Id: rotro2qr.m 2186 2012-07-20 13:45:03Z dmb $
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

% in the comments below, t is the rotation angle, a is the rotation axis
q=zeros(4,1);
t = 1 + r(1) + r(5) + r(9);     % 2(1+cos(t)) = 4(cos(t/2))^2 = 4 q(1)^2
if t>1                       % for rotation angles less than 120 degrees
    s = sqrt(t) * 2;            % 4 cos(t/2) = 2 sin(t)/sin(t/2)
    q(2) = ( r(6) - r(8) ) / s;
    q(3) = ( r(7) - r(3) ) / s;
    q(4) = ( r(2) - r(4) ) / s;
    q(1) = 0.25 * s;            % cos(t/2)
elseif ( r(1) > r(5) && r(1) > r(9) )            % row 1
    s  = sqrt( 1.0 + r(1) - r(5) - r(9) ) * 2;  % 4 a(1) sin(t/2) = 2 a(2) (1-cos(t))/sin(t/2)
    q(2) = 0.25 * s;
    q(3) = (r(2) + r(4) ) / s;
    q(4) = (r(7) + r(3) ) / s;
    q(1) = (r(6) - r(8) ) / s;
elseif ( r(5) > r(9) )                           % row 2
    s  = sqrt( 1.0 + r(5) - r(1) - r(9) ) * 2;
    q(2) = (r(2) + r(4) ) / s;
    q(3) = 0.25 * s;
    q(4) = (r(6) + r(8) ) / s;
    q(1) = (r(7) - r(3) ) / s;
else                                             % row 3
    s  = sqrt( 1.0 + r(9) - r(1) - r(5) ) * 2;
    q(2) = (r(7) + r(3) ) / s;
    q(3) = (r(6) + r(8) ) / s;
    q(4) = 0.25 * s;
    q(1) = (r(2) - r(4) ) / s;
end
f=find(q~=0);
if (q(f(1))<0), q=-q; end           % force positive
if abs(q'*q-1)>1e-8
    error('Input to rotro2qr must be a rotation matrix with det(r)=+1');
end
     