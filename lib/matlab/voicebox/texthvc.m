function h=texthvc(x,y,t,p,q,r)
%TEXTHVC - write text on graph with specified alignment and colour
%
%Usage: (1) texthvc(x,y,'Hello','mlr')          % align centre-baseline colour=red
%       (2) texthvc(x,y,'\alpha','ml',[1 0 0])  % align centre-baseline colour=red
%       (3) texthvc(x,y,'Hello','Mlr')          % x position is normalized to (0,1)
%       (4) texthvc([x k],y,'Hello','mlr')      % add k*axis-width onto the x position
%
% Inputs:  x  x-position of text in graph coordinates (or normalized see below)
%             alternatively x=[x0 k] positions at x0 + k*axis-width
%          y  y-position of text in graph coordinates (or normalized see below)
%             alternatively y=[y0 k] positions at y0 + k*axis-height
%          t  text string to write on graph. Use \alpha for greek or x_1 for subscript
%          p  3-character text string, 'hvc' specifying:
%               horizontal reference point, h: l=left, c or m=middle, r=right
%               vertical reference point,   v: t=top, p=cap, c or m=middle, l=baseline, b=bottom
%               colour,                     c: rgbcmykw
%          q  alterntive colour specification as [r g b] each in range 0 to 1
%
% If the horizontal or vertical reference point is given as a capital
% letter, the corresponding position is normalized to the axis range
% and should be in the range 0 to 1.

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: texthvc.m 8099 2016-06-16 14:48:07Z dmb $
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

if nargin<3
    error('<3 arguments');
elseif nargin>6
    error('>6 arguments');
end

switch nargin
    case 3
        h=text(x,y,char(t));
    case 6
        % for compatibility with old version
        h=text(x,y,char(t),'HorizontalAlignment',p,'VerticalAlignment',q,'Color',r);
    otherwise
        if nargin==5
            r=q; % color specification
        elseif nargin==4
            r=p(3:end);  % colour is at end of spec string
        end
        ix=find(lower(p(1))=='lcmr',1);
        if isempty(ix)
            error('invalid horizontal spec');
        end
        iy=find(lower(p(2))=='tpcmlb',1);
        if isempty(iy)
            error('invalid vertical spec');
        end
        vx={'left', 'center', 'center', 'right'};
        vy={'top', 'cap', 'middle', 'middle', 'baseline', 'bottom'};
        if numel(x)>1
            if strcmpi(get(gca,'xscale'),'log')
                x=exp(log(x(1))+log(get(gca,'xlim'))*[-1; 1]*x(2));
            else
                x=x(1)+get(gca,'xlim')*[-1; 1]*x(2);
            end
        else
            if p(1)==upper(p(1))
                if strcmpi(get(gca,'xscale'),'log')
                    x=exp(log(get(gca,'xlim'))*[1-x; x]);
                else
                    x=get(gca,'xlim')*[1-x; x];
                end
            end
        end
        if numel(y)>1
            if strcmpi(get(gca,'yscale'),'log')
                y=exp(log(y(1))+log(get(gca,'ylim'))*[-1; 1]*y(2));
            else
                y=y(1)+get(gca,'ylim')*[-1; 1]*y(2);
            end
        else
            if p(2)==upper(p(2))
                if strcmpi(get(gca,'yscale'),'log')
                    y=exp(log(get(gca,'ylim'))*[1-y; y])
                else
                    y=get(gca,'ylim')*[1-y; y];
                end
            end
        end
        h=text(x,y,char(t),'HorizontalAlignment',vx{ix},'VerticalAlignment',vy{iy},'Color',r);
end