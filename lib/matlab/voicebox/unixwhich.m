function f=unixwhich(c,e)
%UNIXWHICH Search system path for an executable program [F]=(C,E)
%
% Inputs: C  name of file to search for (excluding extension)
%         E  list of extensions [default = '.com;.exe;.bat' unless C contains '.']
%
% Outputs: F  Full pathname of executable (use FILEPARTS() to split up)
%
% Notes: (1) This is WINDOWS specific and needs to be fixed to work on UNIX systems
%        (2) The search is case insensitive (like most of WINDOWS). 
%        (3) The routine does not cache the directory listings so you
%            should avoid doing the same search many times if you care
%            about speed.
%        (4) To include all files that CMD.EXE will search, set e=winenvar('pathext')
%        (5) As well as their normal full-length name, WINDOWS files and folders have 
%            a short name assigned by the operating system that is 8 characters long
%            (+ 3 more for the extension). These short names are usually hidden from the
%            user and UNIXWHICH, unlike the operating system, will not search for them.

%   Copyright (c) 2005 Mike Brookes,  mike.brookes@ic.ac.uk
%      Version: $Id: unixwhich.m 713 2011-10-16 14:45:43Z dmb $
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
    if any(c=='.')
        e=[];
    else
        e='.com;.exe;.bat';
    end
end
ei=[0 find(e==';') length(e)+1];
[v,pth]=system('path');
while pth(end)<=' '
    pth(end)=[];    % remove junk from end
end
lpth=length(pth);
sc=[0 find(pth==';') lpth+1];
f=[];   % initialize to null string = not found
for i=2:length(sc)
    hd=pth(sc(i-1)+1:sc(i)-1);
    if length(hd)
        [v,fl]=system(['dir /B "',hd,'"']);
        fi=[0 find(fl==10)]; % split into individual file names using LF character
        for j=2:length(fi)
            fn=fl(fi(j-1)+1:fi(j)-1);
            for k=2:length(ei)
                ma=strcmpi(fn,[c e(ei(k-1)+1:ei(k)-1)]);
                if ma
                    f=fullfile(hd,fn);
                    return;
                end
            end
        end
    end
end
