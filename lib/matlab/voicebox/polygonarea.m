function a=polygonarea(p)
%POLYGONAREA Calculate the area of a polygon
%
% Inputs:
%    P(n,2) is the polygon vertices
%
% Outputs:
%    A is teh area of teh polygon
%
% The area is positive if the vertices go anti-clockwise around the polygon.
%

%      Copyright (C) Mike Brookes 2009
%      Version: $Id: polygonarea.m 713 2011-10-16 14:45:43Z dmb $
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
p(end+1,:)=p(1,:);      % append an extra point
a=0.5*sum((p(1:end-1,1)-p(2:end,1)).*(p(1:end-1,2)+p(2:end,2)),1);
if ~nargout
    plot(p(:,1),p(:,2),'b-');
    mnx=[1.05 -0.05;-0.05 1.05]*[min(p);max(p)];
    set(gca,'xlim',mnx(:,1)','ylim',mnx(:,2)');
    title(sprintf('Area = %.2g',a));
end
