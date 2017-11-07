function d=winenvar(n)
%WINENVAR get windows environment variable [D]=(N)
%
% Inputs: N  name of environment variable (e.g. 'temp')
%
% Outputs: D  value of variable or [] is non-existant
%
% Notes: (1) This is WINDOWS specific and needs to be fixed to work on UNIX
%        (2) The search is case insensitive (like most of WINDOWS).
%
% Examples: (1) Open a temporary text file:
%               d=winenar('temp'); fid=fopen(fullfile(d,'temp.txt'),'wt');

%   Copyright (c) 2005 Mike Brookes,  mike.brookes@ic.ac.uk
%      Version: $Id: winenvar.m 713 2011-10-16 14:45:43Z dmb $
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
p=['%',n,'%'];
[s,d]=system(['echo ',p]);
while d(end)<=' ';
    d(end)=[];
end
if strcmp(d,p)
    d=[];
end