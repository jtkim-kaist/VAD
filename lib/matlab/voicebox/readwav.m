function [y,fs,wmode,fidx]=readwav(filename,mode,nmax,nskip)
%READWAV  Read a .WAV format sound file [Y,FS,WMODE,FIDX]=(FILENAME,MODE,NMAX,NSKIP)
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
%             'g'    Scale by the gain written by the "g" option in "writewav" to restore original level
%     Offset: 'y' *	 Correct for offset in <=8 bit PCM data
%             'z'    No offset correction
%   File I/O: 'f'    Do not close file on exit
%             'd'    Look in data directory: voicebox('dir_data')
%   Display;  'h'    Print header information
%             'w'    Plot waveform
%             'W'    Plot spectrogram
%             'a'    play audio (max 10 seconds)
%             'A'    play all audio even if very long
%
%	NMAX     maximum number of samples to read (or -1 for unlimited [default])
%	NSKIP    number of samples to skip from start of file
%               (or -1 to continue from previous read when FIDX is given instead of FILENAME [default])
%
% Output Parameters:
%
%	Y        data matrix of dimension (samples,channels)
%	FS       sample frequency in Hz
%	WMODE    mode string needed for WRITEWAV to recreate the data file
%	FIDX     Information row vector containing the element listed below.
%
%           (1)  file id
%		    (2)  current position in file
%           (3)  dataoff	byte offset in file to start of data
%           (4)  nsamp	number of samples
%           (5)  nchan	number of channels
%           (6)  nbyte	bytes per data value
%           (7)  bits	number of bits of precision
%           (8)  code	Data format: 1=PCM, 2=ADPCM, 3=floating point, 6=A-law, 7=Mu-law
%           (9)  fs	    sample frequency
%           (10) mask   channel mask
%           (11) gain   gain factor in dB
%
%   If no output parameters are specified, header information will be printed.
%
%   For stereo data, y(:,1) is the left channel and y(:,2) the right
%   The mask, if specified, is a bit field giving the channels present in the following order:
%   0=FL, 1=FR, 2=FC, 3=W, 4=BL, 5=BR, 6=FLC, 7=FRC, 8=BC, 9=SL, 10=SR, 11=TC, 12=TFL, 13=TFC, 14=TFR, 15=TBL, 16=TBC, 17=TBR
%   where F=front, L=left, C=centre, W=woofer (low frequency), B=back, LC=left of centre, RC=right of centre, S=side, T=top
%
%   See also WRITEWAV.

%   *** Note on scaling ***
%   If we want to scale signal values in the range +-1 to an integer in the
%   range [-128,127] then we have four plausible choices corresponding to
%   scale factors of (a) 127, (b) 127.5, (c) 128 or (d) 128.5 but each choice
%   has disadvantages.
%   For forward scaling: (c) and (d) cause clipping on inputs of +1.
%   For reverse scaling: (a) and (b) can generate output values < -1.
%   Any of these scalings can be selected via the mode input: (a) 'o', (b) default, (c) 'on', (d) 'n'

%	   Copyright (C) Mike Brookes 1998-2011
%      Version: $Id: readwav.m 713 2011-10-16 14:45:43Z dmb $
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

% Bugs/suggestions:


if nargin<1
    error('Usage: [y,fs,wmode,fidx]=READWAV(filename,mode,nmax,nskip)'); end
if nargin<2
    mode='p';
else
    mode = [mode(:).' 'p'];
end
k=find((mode>='p') & (mode<='s'));
mno=all(mode~='o');                      % scale to input limits not output limits
sc=mode(k(1));
z=128*all(mode~='z');
info=zeros(1,11);
if ischar(filename)
    if any(mode=='d')
        filename=fullfile(voicebox('dir_data'),filename);
    end
    fid=fopen(filename,'rb','l');
    if fid == -1
        fn=[filename,'.wav'];
        fid=fopen(fn,'rb','l');
        if fid ~= -1
            filename=fn;
        end
    end
    if fid == -1
        error('Can''t open %s for input',filename);
    end
    info(1)=fid;
else
    info=filename;
    fid=info(1);
end
getdat= nargout>0 || any(lower(mode)=='w') || any(lower(mode)=='a');
mh=any(mode=='h') || ~getdat;
if ~info(3)
    fseek(fid,8,-1);						% read riff chunk
    header=fread(fid,4,'*char')';
    if ~strcmp(header,'WAVE')
        fclose(fid);
        error('File does not begin with a WAVE chunck');
    end
    if mh
        fprintf('\nWAVE file: %s\n',filename);
    end
    fmtlen=-1;
    datalen=-1;
    instlen=-1;
    factlen=-1;
    riffmt='e';  % default is original wave file format
    while datalen<0						% loop until FMT and DATA chuncks both found
        header=fread(fid,4,'*char');
        len=fread(fid,1,'ulong');
        if mh
            fprintf('  %s chunk: %d bytes\n',header,len);
        end
        if strcmp(header','fmt ')					% ******* found FMT chunk *********
            fmtlen=len;                             % remember the length
            if len>16
                riffmt='x';  % might be WAVEFORMATEX format
            end
            wavfmt=fread(fid,1,'short');			% format: 1=PCM, 6=A-law, 7-Mu-law
            info(8)=wavfmt;
            info(5)=fread(fid,1,'ushort');			% number of channels
            fs=fread(fid,1,'ulong');				% sample rate in Hz
            info(9)=fs;                             % sample rate in Hz
            rate=fread(fid,1,'ulong');				% average bytes per second (ignore)
            align=fread(fid,1,'ushort');			% block alignment in bytes (container size * #channels)
            bps=fread(fid,1,'ushort');			% bits per sample
            info(7)=bps;
            %             info(6)=ceil(info(7)/8);                % round up to a byte count
            info(6)=floor(align/info(5));                       % assume block size/channels = container size
            if info(8)==-2   % wave format extensible
                cb=fread(fid,1,'ushort');			% extra bytes must be >=22
                riffmt='X';                         % WAVEFORMATEXTENSIBLE format
                wfxsamp=fread(fid,1,'ushort');			% samples union
                if wfxsamp>0
                    info(7)=wfxsamp;     % valid bits per sample
                end
                info(10)=fread(fid,1,'ulong');				% channel mask
                wfxguida=fread(fid,1,'ulong');				% GUID
                wfxguidb=fread(fid,2,'ushort');				% GUID
                wfxguidc=fread(fid,8,'uchar');				% GUID
                if wfxguida<65536
                    info(8)=wfxguida; % turn it into normal wav format
                end
                fseek(fid,len-40,0);				    % skip to end of header
            else
                if align>0 && align<(info(6)+4)*info(5)
                    info(6)=ceil(align/info(5));
                end
                fseek(fid,len-16,0);				    % skip to end of header
            end
            if mh
                fmttypes={'?' 'PCM' 'ADPCM' 'IEEE-float' '?' '?' 'A-law' 'µ-law' '?'};
                fprintf('        Format: %d = %s',info(8),fmttypes{1+max(min(info(8),8),0)});
                if wavfmt==-2
                    fprintf(' (%08x-%04x-%04x-%02x%02x-%02x%02x%02x%02x%02x%02x)\n',wfxguida,wfxguidb,wfxguidc);
                else
                    fprintf('\n');
                end
                fprintf('        %d channels at %g kHz sample rate (%d kbytes/s)\n',info(5),fs/1000,rate/1000);
                fprintf('        Mask=%x:',info(10));
                spkpos={'FL' 'FR' 'FC' 'W' 'BL' 'BR' 'FLC' 'FRC' 'BC' 'SL' 'SR' 'TC' 'TFL' 'TFC' 'TFR' 'TBL' 'TBC' 'TBR'};
                for i=1:18
                    if mod(floor(info(10)*pow2(1-i)),2)
                        fprintf([' ' spkpos{i}]);
                    end
                end
                fprintf('\n        %d valid bits of %d per sample (%d byte block size)\n',info(7),bps,align);
            end
        elseif strcmp(header','fact')				% ******* found FACT chunk *********
            factlen=len;
            if len<4
                error('FACT chunk too short');
            end
            nsamp=fread(fid,1,'ulong');				% number of samples
            fseek(fid,len-4,0);                     % skip to end of header
            if mh
                fprintf('        %d samples = %.3g seconds\n',nsamp,nsamp/fs);
            end
        elseif strcmp(header','inst')				% ******* found INST chunk *********
            instlen=len;
            if len<7
                error('INST chunk too short');
            end
            inst=fread(fid,3,'schar');
            info(11)=double(inst(3));                          % gain in dB
            if mh
                fprintf('        Gain = %d dB\n',info(11));
            end
            fseek(fid,len-3,0);                     % skip to end of header
        elseif strcmp(header','data')				% ******* found DATA chunk *********
            if fmtlen<0
                fclose(fid);
                error('File %s does not contain a FMT chunck',filename);
            end
            if factlen>3 && nsamp >0
                info(4)=nsamp;   % take data length from FACT chunk
            else
                info(4) = fix(len/(info(6)*info(5)));  % number of samples
            end
            info(3)=ftell(fid);                    % start of data
            datalen=len;
            if mh
                fprintf('        %d samples x %d channels x %d bytes/samp',info(4:6));
                if prod(info(4:6))~=len
                    fprintf(' + %d padding bytes',len-prod(info(4:6)));
                end
                fprintf(' = %g sec\n',info(4)/fs);
            end
        else							% ******* found unwanted chunk *********
            fseek(fid,len,0);
        end
    end
else
    fs=info(9);
end
if nargin<4 || nskip<0
    nskip=info(2);  % resume at current file position
end

ksamples=info(4)-nskip; % number of samples remaining
if nargin>2
    if nmax>=0
        ksamples=min(nmax,ksamples);
    end
elseif ~getdat
    ksamples=min(5,ksamples); % just read 5 samples so we can print the first few data values
end
if ksamples>0
    info(2)=nskip+ksamples;
    fseek(fid,info(3)+info(6)*info(5)*nskip,-1);
    nsamples=info(5)*ksamples;
    if any(info(8)==3)  % floating point format
        pk=1;  % peak is 1
        switch info(6)
            case 4
                y=fread(fid,nsamples,'float32');
            case 8
                y=fread(fid,nsamples,'float64');
            otherwise
                error('cannot read %d-byte floating point numbers',info(6));
        end
    else
        if ~any(info(8)==[1 6 7])
            sc='r';  % read strange formats in raw integer mode
        end
        pk=pow2(0.5,8*info(6))*(1+(mno/2-all(mode~='n'))/pow2(0.5,info(7)));  % use modes o and n to determine effective peak
        switch info(6)
            case 1
                y=fread(fid,nsamples,'uchar');
                if info(8)==1
                    y=y-z;
                elseif info(8)==6
                    y=pcma2lin(y,213,1);
                    pk=4032+mno*64;
                elseif info(8)==7
                    y=pcmu2lin(y,1);
                    pk=8031+mno*128;
                end
            case 2
                y=fread(fid,nsamples,'short');
            case 3
                y=fread(fid,3*nsamples,'uchar');
                y=reshape(y,3,nsamples);
                y=([1 256 65536]*y-pow2(fix(pow2(y(3,:),-7)),24)).';
            case 4
                y=fread(fid,nsamples,'long');
            otherwise
                error('cannot read %d-byte integers',info(6));
        end
    end
    if sc ~= 'r'
        if sc=='s'
            sf=1/max(abs(y(:)));
        elseif sc=='p'
            sf=1/pk;
        else
            if info(8)==7
                sf=2.03761563/pk;
            else
                sf=2.03033976/pk;
            end
        end
        y=sf*y;
    else                             % mode = 'r' - output raw values
        if info(8)==1
            y=y*pow2(1,info(7)-8*info(6));  % shift to get the bits correct
        end
    end
    if any(mode=='g') && info(11)~=0
        y=y*10^(info(11)/20);   % scale by the gain
    end
    if info(5)>1
        y = reshape(y,info(5),ksamples).';
    end
else
    y=[];
end
if all(mode~='f')
    fclose(fid);
end
if nargout>2  % sort out the mode input for writing this format
    wmode=char([riffmt sc 'z'-z/128]);
    if factlen>0
        wmode=[wmode 'E'];
    end
    if info(6)>1 && info(6)<5
        cszc=' cCL';
        wmode=[wmode cszc(info(6))];
    end
    switch info(8)
        case 1                                    % PCM modes
            if ~mno
                wmode=[wmode 'o'];
            end
            if any(mode=='n')
                wmode=[wmode 'n'];
            end
            wmode=[wmode num2str(info(7))];
        case 3
            if info(7)<=32
                wmode = [wmode 'v'];
            else
                wmode = [wmode 'V'];
            end
        case 6
            wmode = [wmode 'a'];
        case 7
            wmode = [wmode 'u'];
    end
    fidx=info;
end
[ns,nchan]=size(y);
if mh && ns>0
    nsh=min(ns,5);  % print first few samples
    for i=1:nsh
        fprintf('        %d:',i);
        fprintf(' %.3g',y(i,:));
        fprintf('\n');
    end
end

if ns>0.01*fs
    if any(lower(mode)=='a')
        nsh=min(ns,10*fs+ns*any(mode=='A'));
        soundsc(y(1:nsh,1:min(nchan,2)),fs);
    end
    if any(mode=='W')
        spm='pJcbf ';
        if any(mode=='w')
            spm(end)='w';
        end
        clf;
        if nchan>1
            for i=nchan:-1:1
                subplot(nchan,1,i)
                spgrambw(y(:,i),fs,spm);
            end
        else
            spgrambw(y,fs,spm);
        end
    elseif any(mode=='w')
        clf;
        if nchan>1
            for i=nchan:-1:1
                subplot(nchan,1,i)
                plot((1:ns)/fs,y(:,i));
                ylabel(['Chan ' num2str(i)]);
                if i==nchan
                    xlabel('Time (s)');
                end
            end
        else
            plot((1:ns)/fs,y);
            xlabel('Time (s)');
        end

    end
end




