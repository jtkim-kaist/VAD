function [y,fs]=readflac(filename,mode)
%READWAV  Read a .FLAC format sound file [Y,FS]=(FILENAME,MODE)
%
% Input Parameters:
%
%	FILENAME gives the name of the file (with optional .WAV extension) or alternatively
%                 can be the FIDX output from a previous call to READWAV
%	MODE		specifies the following (*=default):
%
%    Scaling: 's'    Auto scale to make data peak = +-1
%             'r'    Raw unscaled data (integer values)
%             'q'    Scaled to make 0dBm0 be unity mean square
%             'p' *	 Scaled to make +-1 equal full scale
%             'o'    Scale to bin centre rather than bin edge (e.g. 127 rather than 127.5 for 8 bit values)
%                     (can be combined with n+p,r,s modes)
%             'n'    Scale to negative peak rather than positive peak (e.g. 128.5 rather than 127.5 for 8 bit values)
%                     (can be combined with o+p,r,s modes)
%
% FLAC (Free lossless audio codec) is a compressed audio file format described here:
% http://flac.sourceforge.net/

%      Copyright (C) Mike Brookes 2008
%      Version: $Id: readflac.m 8483 2016-09-13 15:48:05Z dmb $
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
    mode='p';
else
    mode = [mode(:).' 'p'];
end
dirt=voicebox('dir_temp');
[fnp,fnn,fne]=fileparts(filename);
filetemp=fullfile(dirt,[fnn '.wav']);
doscom=['"' voicebox('flac') '"' ' -d -f -o "' filetemp '" "' filename '"'];
%                     fprintf(1,'Executing: %s\n',doscom);
[doss,dosr]=dos(doscom); % run the program
if doss % test for errors
    error(sprintf('Error running DOS command: %s',doscom));
end
if exist(filetemp)~=2
    error(sprintf('No output file from: %s',doscom));
end
[y,fs]=readwav(filetemp,mode);
doscom=['del /f "' filetemp '"'];
if dos(doscom) % run the program
    error(sprintf('Error running DOS command: %s',doscom));
end
