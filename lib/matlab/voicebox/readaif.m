function [y,fs,wmode,fidx]=readaif(filename,mode,nmax,nskip)
%READAIF  Read a .AIF format sound file [Y,FS,WMODE,FIDX]=(FILENAME,MODE,NMAX,NSKIP)
%
% Input Parameters:
%
%	FILENAME gives the name of the file (with optional .AIF extension) or alternatively
%                 can be the FIDX output from a previous call to READAIF
%	MODE		specifies the following (*=default):
%
%    Scaling: 's'    Auto scale to make data peak = +-1
%             'r'    Raw unscaled data (integer values)
%             'q'    Scaled to make 0dBm0 be unity mean square
%             'p' *	Scaled to make +-1 equal full scale
%             'o'    Scale to bin centre rather than bin edge (e.g. 127 rather than 127.5 for 8 bit values)
%                     (can be combined with n+p,r,s modes)
%             'n'    Scale to negative peak rather than positive peak (e.g. 128.5 rather than 127.5 for 8 bit values)
%                     (can be combined with o+p,r,s modes)
%   File I/O: 'f'    Do not close file on exit
%             'd'    Look in data directory: voicebox('dir_data')
%
%	NMAX     maximum number of samples to read (or -1 for unlimited [default])
%	NSKIP    number of samples to skip from start of file
%               (or -1 to continue from previous read when FIDX is given instead of FILENAME [default])
%
% Output Parameters:
%
%	Y        data matrix of dimension (samples,channels)
%	FS       sample frequency in Hz
%	WMODE    mode string needed to recreate the data file
%	FIDX     Information row vector containing the element listed below.
%
%           (1)  file id
%				(2)  current position in file
%           (3)  dataoff	byte offset in file to start of data
%           (4)  nsamp	number of samples
%           (5)  nchan	number of channels
%           (6)  nbyte	bytes per data value
%           (7)  bits	number of bits of precision
%           (8)  code	Data format: 1=PCM, 2=ADPCM, 6=A-law, 7=Mu-law (always 1 for AIFF)
%           (9)  fs	sample frequency
%           (10) offset
%           (11) block size
%
%   If no output parameters are specified, header information will be printed.
%
%   For stereo data, y(:,1) is the left channel and y(:,2) the right

%  future enhancements: (1) handle LIST and CAT files as well
%                       (2) deal with offset and blocksize properly
%                       (3) handle other chunk types sensibly

%	   Copyright (C) Mike Brookes 1998
%      Version: $Id: readaif.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<1 error('Usage: [y,fs,wmode,fidx]=READAIFF(filename,mode,nmax,nskip)'); end

if nargin<2 mode='p';
else mode = [mode(:).' 'p'];
end
k=find((mode>='p') & (mode<='s'));
mno=all(mode~='o');                      % scale to input limits not output limits
sc=mode(k(1)); 


info=zeros(1,11);
if ischar(filename)
    if any(mode=='d')
        filename=fullfile(voicebox('dir_data'),filename);
    end
    fid=fopen(filename,'rb','b');
    if fid == -1
        fn=[filename,'.aif'];
        fid=fopen(fn,'rb','b'); 
        if fid ~= -1 filename=fn; end
    end
    if fid == -1 error(sprintf('Can''t open %s for input',filename)); end
    info(1)=fid;
else
    info=filename;
    fid=info(1);
end

if ~info(3)
    fseek(fid,0,-1);						% read file header
    header=fread(fid,4,'uchar');
    if header' ~= 'FORM' fclose(fid); error(sprintf('File %s does not begin with a FORM group ID',filename)); end
    filen=fread(fid,1,'ulong')-4;
    header=fread(fid,4,'uchar');
    if header' ~= 'AIFF' fclose(fid); error(sprintf('File %s does not begin with a AIFF type ID',filename)); end
    
    
    fmt=0;
    data=0;
    while filen>=8						% loop to read all chunks
        header=fread(fid,4,'char');
        len=fread(fid,1,'ulong');
        lenx=len+rem(len,2);    % round up to an even number
        filen=filen-lenx-8;
        %fprintf(1,'%s chunk, %d bytes, %d bytes remaining\n',char(header'),len,filen);
        if header' == 'COMM'					% ******* found COMM chunk *********
            fmt=1;
            info(5)=fread(fid,1,'ushort');			%number of channels
            info(4)=fread(fid,1,'ulong');			%number of samples
            info(7)=fread(fid,1,'ushort');				% bits per sample
            info(6)=ceil(info(7)/8);            % bytes per sample
            info(8)=1;
            ieeex=fread(fid,5,'ushort');               % read a 10-byte extended ieee format number
            fsign=ieeex(1)>32767;
            ieeex(1)=ieeex(1)-32768*fsign;
            if ~ieeex info(9)=0;
            elseif ieeex(1)==32767 info(9)=nan;
            else
                info(9)=(1-2*fsign)*sum(pow2(ieeex(2:5)',ieeex(1)-(16398:16:16446)));
            end
            fseek(fid,lenx-18,0);				% skip to end of header
        elseif header' == 'SSND'				% ******* found DATA chunk *********
            info(10)=fread(fid,1,'ulong');			%number of channels
            info(11)=fread(fid,1,'ulong');			%number of channels
            info(3)=ftell(fid);
            data=1;
            fseek(fid,lenx-8,0);				% skip to end of chunk
        else							% ******* found unwanted chunk *********
            fseek(fid,lenx,0);
        end
    end
    if filen~=0
        fprintf(2,'READAIF warning: Inconsistent file length - %d extra bytes\n',filen);
    end
end
fs=info(9);
if ~fmt | ~data fclose(fid); error(sprintf('File %s does not contain COMM and SSND chunks',filename)); end


if nargin<4 nskip=info(2);
elseif nskip<0 nskip=info(2);
end

ksamples=info(4)-nskip;
if nargin>2
    if nmax>=0
        ksamples=min(nmax,ksamples);
    end
elseif ~nargout
    ksamples=min(5,ksamples);
end
if ksamples>0
    info(2)=nskip+ksamples;
    pk=pow2(0.5,8*info(6))*(1+(mno/2-all(mode~='n'))/pow2(0.5,info(7)));  % use modes o and n to determine effective peak
    fseek(fid,info(3)+info(6)*info(5)*nskip,-1);
    nsamples=info(5)*ksamples;
    if info(6)<3
        if info(6)<2
            y=fread(fid,nsamples,'schar');
        else
            y=fread(fid,nsamples,'short');
        end
    else
        if info(6)<4
            y=fread(fid,3*nsamples,'uchar');
            y=reshape(y,3,nsamples);
            y=[1 256 65536]*y-pow2(fix(pow2(y(3,:),-7)),24);
        else
            y=fread(fid,nsamples,'long');
        end
    end
    if sc ~= 'r'
        if sc=='s' sf=1/max(max(abs(y)),1);
        elseif sc=='p' sf=1/pk;
        else   sf=2.03033976/pk;
        end
        y=sf*y;
    else
        y=y*pow2(1,info(7)-8*info(6)); % shift to get the data into the LSB end
    end
    
    if info(5)>1 y = reshape(y,info(5),ksamples).'; end
else
    y=[];
end

if all(mode~='f') fclose(fid); end

if nargout>2
    if ~mno wmode=[wmode 'o']; end
    if any(mode=='n') wmode=[wmode 'n']; end
    wmode=[sc num2str(info(7))];
    fidx=info;
elseif ~nargout
    fprintf(1,'\n%d Hz sample rate\n%d channel x %d samples x %d bit = %.3g seconds\n',info([9 5 4]), info(7),info(4)/info(9));
end



