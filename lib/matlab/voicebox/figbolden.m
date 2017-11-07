function figbolden(pos,pv,m)
%FIGBOLDEN embolden, resize and recolour the current figure =(POS,PV,M)
% 
% Inputs: pos = [xmin ymin width height] gives the lower left corner position and the window size in pixels
%               [width height] leaves the lower left corner alone
%               [width] has a standard aspect ratio of 4:3
%               [-width/height] leaves the area unchanged but fixes the aspect ratio
%         pv is a cell name containing attribute-value pairs.
%            default = {'FontName' 'Arial'; 'FontSize' 16; 'LineWidth' 2; 'MarkerSize' 8}
%            Note that the "listfonts" MATLAB command will list the available fonts
%         m is a mode string:
%                'l' lists the changes made
%                'd' append default pv settings (e.g. use [] for second argument)
%                'c' change default colours to improve contrast on a white background
%                    g->[0,0.7,0],c->[0,0.7,0.7],y->[0.83,0.83,0]
%                'x' suppresses all changes
%
% Bug: gives an error message if log axes have been used

%      Copyright (C) Mike Brookes 2003
%      Version: $Id: figbolden.m 9186 2016-12-13 08:51:37Z dmb $
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

ps={'Title' 'XLabel' 'YLabel' 'Children'};
if nargin<3
    m='';
end
if nargin<2 || any(m=='d')
    pv={'FontName' 'Arial'; 'FontSize' 16; 'LineWidth' 2; 'MarkerSize' 8};
end
pp={'Symbol';'Wingdings'};      % protected fonts
if nargin<1
    pos=[];
end
mlist=any(m=='l');  % list changes
mnotx=~any(m=='x'); % do changes
scsz=get(0,'screensize');
if length(pos)
    po=get(gcf,'position');
    if length(pos)>2            % position is specified
        po(1:2)=pos(1:2);
        pos(1:2)=[];      % remove xmin,ymin
    end
    if length(pos)>1
        po(3:4)=pos(1:2);
    else
        if pos(1)>0
            po(3:4)=[1 0.75]*pos(1);
        else
            po(3:4)=[-pos(1) 1]*sqrt(-po(3)*po(4)/pos(1)); % preserve area
        end
    end
    set(gcf,'position',po);
end
if any(m=='c')
    cc='gcy';
    cv=[0,0.7,0;0,0.7,0.7;0.83,0.83,0];
    for ic=1:3
        hlist=findobj(gcf,'color',cc(ic));
        if length(hlist)
            if mnotx
                set(hlist,'color',cv(ic,:));
            end
            if mlist
            for ih=1:length(hlist)
                 fprintf(['change %f Color: ''%c'' -> [%g %g %g]\n'],hlist(ih),cc(ic),cv(ic,:));
            end
            end
        end
    end
end
hlist=get(gcf,'children');
while length(hlist)
    pl=get(hlist(1));
    %fprintf('list length = %d, handle = %f\n',length(hlist),hlist(1));
    for i=1:size(pv,1)
        if isfield(pl,pv{i,1})
            if i>1 || all(~strcmpi(get(hlist(1),pv{i,1}),pp))
                pval=get(hlist(1),pv{i,1});
                if ~all(size(pval)==size(pv{i,2})) || ~all(pval(:) == pv{i,2}(:))
                    if mnotx
                        set(hlist(1),pv{i,1},pv{i,2})
                    end
                    if mlist
                        if ischar(pval)
                            fch='s';
                        else
                            fch='g';
                        end
                        fprintf(['change %f %s: %' fch ' -> %' fch '\n'],hlist(1),pv{i,1},pval,pv{i,2});
                    end
                end
            end
        end
    end
    for i=1:length(ps)
        if isfield(pl,ps{i})
            hlist=[hlist; get(hlist(1),ps{i})];
            %fprintf('add %f:%s\n',hlist(1),ps{i});
        end
    end
    hlist(1)=[];
end

