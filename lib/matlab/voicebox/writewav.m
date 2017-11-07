function fidx=writewav(d,fs,filename,mode,nskip,mask)
%WRITEWAV Creates .WAV format sound files FIDX=(D,FS,FILENAME,MODE,NSKIP,MASK)
%
%   The input arguments for WRITEWAV are as follows:
%
%       D           The sampled data to save
%       FS          The rate at which the data was sampled
%       FILENAME    A string containing the name of the .WAV file to create or
%                        alternatively the FIDX output from a previous writewav call
%       MODE        String containing any reasonable mixture of flags below (*=default):
%       NSKIP       Number of samples to skip before writing or -1[default] to continue from previous write
%                   Only valid if FIDX is specified for FILENAME
%       MASK        specifies the speaker positions included as a bit mask (see readwav)
%
% MODE flags (*=default):
%  Precision: 'a'    for 8-bit A-law PCM
%             'u'    for 8-bit mu-law PCM
%            '16' *	 for 16 bit PCM data
%             '8'    for 8 bit PCM data
%             ...    any number in the range 2 to 32 for PCM
%             'v'    32-bit floating point
%             'V'    64-bit floating point
%             'c'    embed in 16 bits
%             'C'    embed in 24 bits
%             'L'    embed in 32 bits
%	  Dither: 'w'    White triangular dither of amplitude +-1 LSB (PCM modes only)
%             'h'    High pass dither (filtered by 1-1/z) (PCM modes only)
%             'l'    Low pass dither (filtered by 1+1/z) (PCM modes only)
%    Scaling: 's' *  Auto scale to make data peak = +-1
%             'r'    Raw unscaled data (integer values)
%             'q'    Scaled to make unity mean square correspond to 0dBm according to G.711
%             'p'  	 Scaled to make +-1 equal full scale
%             'o'    Scale to bin centre rather than bin edge (e.g. 127 rather than 127.5 for 8 bit values)
%                     (can be combined with n+p,r,s modes)
%             'n'    Scale to negative peak rather than positive peak (e.g. 128.5 rather than 127.5 for 8 bit values)
%                     (can be combined with o+p,r,s modes)
%             'g'    Include a gain factor so that "readwav" will restore the correct level
%     Offset: 'y' *	 Correct for offset in <=8 bit PCM data
%             'z'    Do not apply offset correction
%     Format: 'x'    use WAVEFORMATEX format (default for non PCM)
%             'X'    use WAVEFORMATEXTENSIBLE (default if MASK input is given)
%             'e'    use original WAVEFORMAT (default for PCM)
%             'E'    include a 'fact' chunk (default for non-PCM)
%   File I/O: 'f'    Do not close file on exit
%             'd'    Look in data directory: voicebox('dir_data')
%
%
% Output Parameter:
%
%	FIDX     Information row vector containing the element listed below.
%
%           (1)  file id
%			(2)  current position in file (in samples, 0=start of file)
%           (3)  dataoff	length of file header in bytes
%           (4)  nsamp	number of samples
%           (5)  nchan	number of channels
%           (6)  nbyte	bytes per data value
%           (7)  bits	number of bits of precision
%           (8)  code	Data format: 1=PCM, 2=ADPCM, 6=A-law, 7=Mu-law
%           (9)  fs	sample frequency
%           (10) dither state variable
%           (11) gain in dB (in INST chunk)
%
%   Note: WRITEWAV will create an 16-bit PCM, auto-scaled wave file by default.
%   For stereo data, d(:,1) is the left channel and d(:,2) the right
%
%   See also READWAV

%   *** Note on scaling ***
%   If we want to scale signal values in the range +-1 to an integer in the
%   range [-128,127] then we have four plausible choices corresponding to
%   scale factors of (a) 127, (b) 127.5, (c) 128 or (d) 128.5 but each choice
%   has disadvantages.
%   For forward scaling: (c) and (d) cause clipping on inputs of +1.
%   For reverse scaling: (a) and (b) can generate output values < -1.
%   Any of these scalings can be selected via the mode input: (a) 'o', (b) default, (c) 'on', (d) 'n'

%	   Copyright (C) Mike Brookes 1998-2011
%      Version: $Id: writewav.m 713 2011-10-16 14:45:43Z dmb $
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

% Acknowledgements
%   Thanks to Hugh Barnes for sorting out seek problems with MATLAB 6.5

% Bugs/suggestions
%  Save the following factors in FIDX: (a) scale factor, (b) offset (c) low/high clip limits
%       (d) dither position  (e) rounding position


if nargin<3
    error('Usage: WRITEWAV(data,fs,filename,mode,nskip)');
end
if nargin<6
    mask=0;
end
info=zeros(1,11);
info(9)=fs;
if nargin<4
    mode='s';
else
    mode = [mode(:).' 's'];  % ensure that there is always a scaling option specified
end
info(8)=1;     % default mode is PCM
mno=all(mode~='o');                      % scale to input limits not output limits
k=find((mode>='0') & (mode<='9'),1);
if k,
    info(7)=sscanf(mode(k:end),'%d');  % valid bits per data value
else
    info(7)=16;
end
if any(mode=='c')
    info(6)=2;       % bytes per data value = 2
elseif any(mode=='C')
    info(6)=3;       % bytes per data value = 3
elseif any(mode=='L')
    info(6)=4;       % bytes per data value = 4
else
    info(6)=ceil(info(7)/8);       % bytes per data value
end
lo=-pow2(0.5,info(7));
hi=-1-lo;
pk=pow2(0.5,8*info(6))*(1-(mno/2-all(mode~='n'))/lo);  % use modes o and n to determine effective peak
% should perhaps have another variable besides info(7) to control dither position (or set info(7) later)
% for A and mu law the dither position is not the same as the number of bits.
if any(mode=='a')
    info(8)=6;
    pk=4032+mno*64;
    info(7)=8;  % Some sources say this should be listed as 16 valid bits
    info(6)=1;
elseif any(mode=='u')
    info(8)=7;
    pk=8031+mno*128;
    info(7)=8;  % Some sources say this should be listed as 16 valid bits
    info(6)=1;
elseif any(mode=='v')
    pk=1;
    mode(end)='r';  % default scaling is 'r'
    info(6)=4;  % bytes
    info(7)=32; % bits
    info(8)=3;  % WAVE type
elseif any(mode=='V')
    pk=1;
    mode(end)='r';   % default scaling is 'r'
    info(6)=8; % bytes
    info(7)=64; % bits
    info(8)=3; % WAVE type
end			% is this pk value correct ?
sc=mode(find((mode>='p') & (mode<='s'),1)); % find the first scaling option (always exists)
z=128*all(mode~='z');
if any(mode=='w')
    di='w';                       % select dither mode
elseif any(mode=='h')
    di='h';
elseif any(mode=='l')
    di='l';
else
    di='n';
end

% Now sort out which wave format to use
if any(mode=='e')
    wavtype=1;
elseif any(mode=='x')
    wavtype=2;
elseif any(mode=='X') || nargin>=6
    wavtype=3;
else
    wavtype=2-(info(8)==1);
end
wavfmt=info(8)*(wavtype<3)+(pow2(16)-2)*(wavtype==3);
fmtlen=[16 18 40]; % length of format chunk
factlen=12*(any(mode=='E') || info(8)~=1);
instlen=16*any(mode=='g');  % length of INST chunk (force to be even since some readers do not like odd lengths)
wavlen=[36 38 60]+factlen+instlen; % length of entire WAVE chunk except for the data (not including 8 byte RIFF header)


[n,nc]=size(d);
if n==1
    n=nc;
    nc=1;
else
    d = d.';
end;
if nc>32
    error('WRITEWAV: attempt to write a sound file with >32 channels');
end
nc=max(nc,1);
ncy=nc*info(6);                     % bytes per sample time
nyd=n*ncy;                          % bytes to write

if ischar(filename)
    if any(mode=='d')
        filename=fullfile(voicebox('dir_data'),filename);
    end
    ny=nyd;
    if isempty(findstr(filename,'.'))
        filename=[filename,'.wav'];
    end
    fid=fopen(filename,'wb+','l');
    if fid == -1
        error('Can''t open %s for output',filename);
    end
    info(1)=fid;
    fwrite(fid,'RIFF','uchar');  % main RIFF header
    fwrite(fid,wavlen(wavtype)+2*ceil(ny/2),'ulong');  %
    fwrite(fid,'WAVEfmt ','uchar');   % write "WAVE" ID and "fmt" chunk
    fwrite(fid,[fmtlen(wavtype) 0 wavfmt nc],'ushort'); % chunk size, format code, number of channels
    fwrite(fid,[fs fs*ncy],'ulong');        % sample rate, bytes per sec
    switch wavtype
        case 1
            fwrite(fid,[ncy info(7)],'ushort');     % block size, bits-per-sample
        case 2
            fwrite(fid,[ncy info(7)],'ushort');     % block size, bits-per-sample
            fwrite(fid,0,'ushort');     % size of the extension=0
        case 3
            fwrite(fid,[ncy 8*info(6)],'ushort');     % block size, bits-per-sample (aways a multiple of 8)
            fwrite(fid,[22 info(7)],'ushort');     % size of the extension=22, valid bits
            fwrite(fid,[mask info(8)],'ulong');        % speaker position mask, encoding format
            fwrite(fid,[0 16 128 43520 14336 29083],'ushort');				% GUID
    end
    if factlen
        fwrite(fid,'fact','uchar');   % fact chunk header
        fwrite(fid,[4 n],'ulong');       % length in bytes + number of samples
    end
    if instlen
        fwrite(fid,'inst','uchar');   % fact chunk header
        fwrite(fid,instlen-8,'ulong');       % length in bytes
        fwrite(fid,zeros(1,instlen-8),'uchar');   % inst data (zero for now)
    end
    fwrite(fid,'data','uchar');   % data header
    fwrite(fid,ny,'ulong');       % data length in bytes
    nskip=0;                        % over-ride any nskip argument
    info(3)=8+wavlen(wavtype);      % length of all header information
    info(4)=n;                      % number of samples (per channel)
    info(2)=n;                      % current file position (in samples)
    info(10)=rand(1);                       % seed for dither generation
else
    info=filename;
    fid=info(1);
    fseek(fid,0,1); % go to end of file
    if nargin<5 || nskip<0
        nskip=info(2);                      % use previous high water mark
    end
    info(2)=n+nskip;                      	% file position following this write operation (in samples)
    ny=nyd+nskip*ncy;                       % file position following this write operation (in bytes following header)
    if n && (info(2)>info(4))               % update high water mark
        if ~info(4)                       	% if no data written previously
            fseek(fid,22,-1); fwrite(fid,nc,'ushort');   % update number of channels
            fseek(fid,28,-1); fwrite(fid,fs*ncy,'ulong'); % update bytes/second
            fwrite(fid,ncy,'ushort'); % update bytes/sample
        end
        fseek(fid,4,-1); fwrite(fid,wavlen(wavtype)+2*ceil(ny/2),'ulong'); % update RIFF length
        if factlen
            fseek(fid,wavlen(wavtype)-4-instlen,-1); fwrite(fid,n,'ulong');  % update FACT number of samples
        end
        fseek(fid,4+wavlen(wavtype),-1); fwrite(fid,ny,'ulong');  % update DATA length
        info(4)=info(2);
    end
end
info(5)=nc;

if n

    if sc~='r'                  % 'r' = no scaling
        if sc=='s'              % 's' = scale to peak signal
            pd=max(abs(d(:)));
            pd=pd+(pd==0);      % scale to 1 if data is all zero
        elseif sc=='p'          % 'p' = scale to +-1 = full scale
            pd=1;
        else                    % 'q' = scale to 0dBm
            if info(8)==7       % mu-law
                pd=2.03761563;
            else                % A-law or anything else
                pd=2.03033976;
            end
        end
        if instlen
            info(11)=min(max(ceil(20*log10(pd)),-128),127);
            d=pk*10^(-0.05*info(11))*d;    
            if fseek(fid,0,-1)  % MATLAB V6.5 fails if this is omitted
                error('Cannot rewind file');
            end
            if fseek(fid,info(3)-instlen+2,-1);
                error('Cannot seek to INST chunk gain byte');
            end
            fwrite(fid,info(11),'schar');   % write the INST gain in dB
        else
            d=pk/pd*d;
        end
    end
    if fseek(fid,0,-1)  % MATLAB V6.5 fails if this is omitted
        error('Cannot rewind file');
    end
    if fseek(fid,info(3)+nskip*nc*info(6),-1)
        error('Cannot seek to byte %d in output file',info(3)+nskip*nc*info(6));
    end
    if info(8)==3 % floating point
        if info(6)==4
            fwrite(fid,d,'float32');
        else
            fwrite(fid,d,'float64');
        end
    else                        % integer data
        if info(8)<6            % PCM
            if di=='n'
                d=round(d);
            else
                [d,info(10)]=ditherq(d,di,info(10));
            end
            d=min(max(d,lo),hi)*pow2(1,8*info(6)-info(7));       % clip data and shift to most significant bits
        else                    % mu or A law
            z=0;
            if info(8) < 7
                d=lin2pcma(d,213,1);
            else
                d=lin2pcmu(d,1);
            end
        end
        if info(6)<3
            if info(6)<2
                fwrite(fid,d+z,'uchar');
            else
                fwrite(fid,d,'short');
            end
        else
            if info(6)<4
                d=d(:)';
                d2=floor(d/65536);
                d=d-65536*d2;
                fwrite(fid,[rem(d,256); floor(d/256); d2+256*(d2<0)],'uchar');
            else
                fwrite(fid,d,'long');
            end
        end
        if rem(ny,2) % pad to an even number of bytes
            fwrite(fid,0,'uchar');
        end
    end
end
if all(mode~='f')
    fclose(fid);
end
if nargout
    fidx=info;
end
