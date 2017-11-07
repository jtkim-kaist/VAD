function c=cblabel(l,h)
% h is the handle of the colorbar, axis or figure
%CBLABEL add a label to a colorbar c=(l,h)
%
% Inputs:
%
%     L        Label string for colorbar
%     H        Handle of the colorbar, axis or figure [default = current figure]
%
% Outputs:
%
%     C        Handle of the colorbar

%      Copyright (C) Mike Brookes 2000-2009
%      Version: $Id: cblabel.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<2
    h=gcf;
end
switch get(h,'Type')
    case 'axes'
        if get(h,'Tag')=='Colorbar'
            c=h;
        else
            while get(h,'Type')~='figure'
                h=get(h,'Parent');      % find parent figure
                if h==0
                    error('cannot find parent figure');
                end
            end
            c=findobj(h,'tag','Colorbar');
            if isempty(c)
                error('There is no colour bar on this figure')
            end
            % we could look for the nearest colorbar to the selected axes
            c=c(1);      % for now use the most recently added colorbar
        end
    case 'figure'
        c=findobj(h,'tag','Colorbar');
        if isempty(c)
            error('There is no colour bar on this figure')
        end
        c=c(1);      % use the most recently added colorbar
    otherwise
        error('h argument must be colorbar, axis or figure handle');
end
set(get(c,'ylabel'),'string',l);