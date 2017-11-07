function [y,fs,h]=readcnx(filename,mode)
%READCNX  Read a .CNX format sound file [Y,FS,H]=(FILENAME)
%
% Inputs:
%           filename : character string containing filename (default extension .cnx)
%           mode:    't' = trim leading and trailing silences
%                    'h' = read header only
%             'd'    Look in data directory: voicebox('dir_data')
%
% Outputs:
%           y  : column vector containing waveform
%           fs : sample frequency
%           h  : parameter array:
%                h(1) = number of samples in file
%                h(2) = status: 0=good, 1=bad
%                h(3) = start sample number
%                h(4) = ending sample number
%                h(5) = speaker identification number
%                h(6) = speaker age group
%                h(7) = speaker sex: 0=male, 1=female
%                h(8) = ascii character
%                h(9) = repetition number
%
% This is the format of the BT Connex-S1 alphabet database
% Note: the decoding is not particularly robust and assumes
%       that all headers contain the same sequence of fields

%	   Copyright (C) Mike Brookes 1998
%      Version: $Id: readcnx.m 713 2011-10-16 14:45:43Z dmb $
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

ix=[17 71; 18 0; 19 0; 10 0; 12 0; 13 77; 15 -1; 16 0; ];

if nargin<2  mode='0'; end 
if any(mode=='d')
    filename=fullfile(voicebox('dir_data'),filename);
end
fid=fopen(filename,'rb','l');
if fid == -1
    fn=[filename,'.cnx'];
    fid=fopen(fn,'rb','l');
    if fid ~= -1 filename=fn; end
end
if fid == -1 error(sprintf('Can''t open %s for input',filename)); end
[hdr,n]=fread(fid,512,'uchar');
if n ~= 512 error(sprintf('Can''t read header from connex file %s',filename)); end
del=find(hdr(5:end)=='|')';
fs=sscanf(char(hdr(17:del(1)+3)),'%f');
h=zeros(size(ix,1),1);
for i=1:length(h)
    e=find(hdr(del(ix(i)-1)+5:del(ix(i))+3)=='=');
    if ix(i,2)
        h(i)=hdr(del(ix(i)-1)+e+5);
        if ix(i,2)>0
            h(i)=1-(h(i)==ix(i,2));
        end
    else
        h(i)=sscanf(char(hdr(del(ix(i)-1)+e+5:del(ix(i))+3)),'%d');
    end
end
if any(mode =='h')
    y=[];
elseif any(mode =='t')
    fseek(fid,2*h(2),0);
    [y,n]=fread(fid,h(3)-h(2)+1,'short');
    if n ~= h(3)-h(2)+1 error(sprintf('Error reading data from connex file %s',filename)); end
else
    y=fread(fid,inf,'short');
end
fseek(fid,0,1);
h=[ftell(fid)/2-256; h];
fclose(fid);
