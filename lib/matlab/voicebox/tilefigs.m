function tilefigs(pos)
%TILEFIGS tile current figures
% Inputs: pos(1,4) Virtual screen: [xmin ymin width height] either in
%                  pixels >=1 or normalized in the range 0 to 1
%
% Possible future options:
%   (1) define screen region to place tiles
%   (2) include/exclude window titles (except optionally first row)
%   (3) preserve sizes and/or aspect ratios
%   (4) place disparate sizes in the most efficient way
%   (5) cascade instead of tile
%   (6) superimpose instead of tile
%   (7) cope with dual screens
%   (8) should check that the units have not been changed from "pixels"
%   (9) cope with virtual screen that has a different aspect ratio

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: tilefigs.m 8587 2016-09-23 15:53:27Z dmb $
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

starth=35; % height for start bar
tith=78; % height for window title
winb=8; % width of window border
figl=get(0,'Children');
nf=length(figl);
if isnumeric(figl)  % old versions of MATLAB have numeric figure handles
    figl=sort(figl);
else
    fign=zeros(1,nf);
    for i=1:nf
        fign(i)=figl(i).Number;
    end
    [fign,figj]=sort(fign);
    figl=figl(figj);
end
scr=get(0,'Screensize'); % [xmin xmax width height]
if nargin==0 || isempty(pos)
    scr=scr+starth*[0 1 0 -1]; % leave space for the start bar at the bottom
else
    if all(pos<2) % pos uses normalized units
        pos=[round(scr(1:2)+pos(1:2).*scr(3:4)) round(pos(3:4).*scr(3:4))]; % convert to pixels
    end
    scr=[max(pos(1:2),1) min(pos(3:4),scr(3:4)+scr(1:2)-max(pos(1:2),1))]; % clip to actual screen
end
nfh=1:nf; % possible number of columns
nfv=ceil(nf./nfh); % corresponding number of rows
asp=(floor(scr(3)./nfh)-2*winb)./(floor(scr(4)./nfv)-2*winb-tith); % corresponding aspect ratios
[aa,ia]=min(abs(asp-4/3)); % target aspect ratio is 4/3
nfv=nfv(ia);
nfh=ceil(nf/nfv);
hpix=floor(scr(3)/nfh); % horizontal pixels per graph available incl border
vpix=floor(scr(4)/nfv); % vertical pixels per graph available incl border
% fprintf('scr=[%d %d %d %d], [hpix vpix]=[%d %d]\n',scr,hpix,vpix);
for i=nf:-1:1
    figure(figl(i));
    row=1+floor((i-1)/nfh); % row goes nvf:-1:1
    col=i-(row-1)*nfh;      % within each row col goes nfh:-1:1 except for last row
% fprintf('Fig %d=(%d,%d): [xmin ymin width height] = [%d %d %d %d]\n',i,row,col,[hpix*(col-1)+1+winb scr(4)-vpix*row+1+winb hpix-2*winb vpix-2*winb-tith]);
    set(figl(i),'position',[hpix*(col-1)+1+winb scr(4)+scr(2)-vpix*row+winb hpix-2*winb vpix-2*winb-tith]);
end