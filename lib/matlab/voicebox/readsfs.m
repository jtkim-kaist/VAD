function [y,fs,hd,ffx]=readsfs(ff,ty,sub,mode,nmax,nskip,xpath)
%READSFS  Read a .SFS format sound file [Y,FS,HD,FFX]=(FF,TY,SUB,MODE,NMAX,NSKIP,XPATH)
%
% Usage:  [s,fs]=readsfs(filename,1); % read the last speech item in a file
%
% Input Parameters:
%
%  FF gives the name of the file or alternatively
%                 can be the ffx output from a previous call to readsfs
%  TY gives the type of data item required e.g.:
%              0 Main header, 1 Speech data, 2 Laryngograph, 5 Annotation
%  SUB specifies the instance of type TY in the file: 0 for first (default), -1 for last or else
%      it can specify the start of the processing history field as a string (e.g. 'hqtx')
%  MODE		specifies the following (*=default):
%
%           File I/O:    'f'    Do not close file on exit
%                        'd'    Look in data directory: voicebox('dir_data')
%           Int Format:  'i'  Force integer data to be at least 16 bits
%                              (some sfs files have a header error which falsely indicates 8-bit data)
%           Create item: 'c' Create item if necessary
%           Errors:      'x' Ignore errors
%
%  NMAX     maximum number of samples to read (or -1 for unlimited [default])
%  NSKIP    number of samples or frames to skip from start of file
%               (or -1 to continue from previous read when FFX is given instead of a filename [default])
%  XPATH    (used with 'c' option) gives the full name of the program needed to generate the data or
%           the directory containing it.
%
% Output Parameters:
%
%  Y        data matrix or cell matrix whose format depends on TY:
%        TY=0: empty
%			TY=5: cell array {nf,3} = {position length annotation}
%        TY=1,2: column vector containing data
%        TY=11: data array with one row per frame
%  FS       sample frequency in Hz
%  HD     cell matrix whose format depends on TY:
%        TY=0: cell{14,1}
%              {1} row vector
%                {1}(1) = serial_date (see DATENUM() for format)
%                {1}(2) = file_number
%                {1}(3) = machine_type
%              {2} = File type (= 'UC2')
%              {3} = username of creator
%              {4} = site of creator
%              {5} = source
%              {6} = database
%              {7} = speaker name
%              {8} = session code
%              {9} = session date (as a string)
%             {10} = name of token
%             {11} = token repetition code
%             {12} = recording conditions
%             {13} = archiving details
%             {14} = general comments
%        TY>0: cell{4,1}
%              {1} = (1,14) array:
%                 {1}(1)  = processdate (see DATENUM() for format)
%                 {1}(2)  = datatype: 1=speech, 2=lx, 3=tx cycle lengths, 4=fx freq
%                                     5=annotations, 6=phonetic, 7=synthesiser, 8=words
%                                     9=grey-level, 10=voicing, 11=energy, 12=formants
%                                     13=energy, 14=lpc, 15=markov, 16=acoustic, 17=?,
%                                     18=geometry, 19=aerodynamics, 20=articulatory
%                                     21=source, 22=physiological, 23=rational filter
%                                     24=poles/zeros, 25=glottal flow, 26=excitation model
%                                     27=nose, 28=calibration, 29=area
%                 {1}(3)  = subtype
%                 {1}(4)  = floating: 1=float, 0=int, -1=structure+
%                 {1}(5)  = datasize in bytes
%                 {1}(6)  = framesize in units of datasize
%                 {1}(7)  = numframes
%                 {1}(8)  = length in bytes of data
%                 {1}(9)  = frameduration=1/sample_rate
%                 {1}(10) = datapresent: 0=deleted, 1=present, 2=link
%                 {1}(11) = time offset
%                 {1}(12) = windowsize
%                 {1}(13) = overlap
%                 {1}(14) = lxsync
%              {2} = processing history
%              {3} = parameter field
%              {4} = comment
%
%  FFX     cell array containing:
%              {1} = filename
%              {2} = (1,4) = [fid byte_order item_row values_read]
%              {3} = (nitem,5) = one row per item [type subtype length position byteorder]
%                                with deleted items having a negative type and zero length
%              {4} = {nitem,3} = cell: one row per item {processing parameters comment} text strings
%
%
% The SFS (Speech Filing System) is a package mostly written by Mark Huckvale
% and is available for UNIX and PC systems from
% http://www.phon.ucl.ac.uk/resource/sfs/

% Features yet to be implemented:
%
%   (1) If no output parameters are specified, header information will be printed.
%   (2) following link items
%   (3) MODE options:
%                Scaling: 's'    Auto scale to make data peak = +-1
%                         'r'    Raw unscaled data (integer values)
%                         'q'    Scaled to make 0dBm0 be unity mean square
%                         'p' *	Scaled to make +-1 equal full scale
%                Errors   'r'    Return if file is non-existant


%	   Copyright (C) Mike Brookes 1998
%      Version: $Id: readsfs.m 6803 2015-09-12 09:31:44Z dmb $
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

if nargin<7
    xpath=voicebox('sfsbin');       % path for sfs programs
end
EXESUF=voicebox('sfssuffix');                                              % suffix for executable O/S dependent
if nargin<4 || ~numel(mode)
    mode='p';
else
    mode = [mode(:).' 'p'];
end
if nargout==0
    if nargin<2 || ~numel(ty) || ty<=0
        [yy,fs,hd,ffx]=readsfs(ff,0,0,mode);
        fprintf(1,'File: %s\n',ffx{1});
        fprintf(1,'Database: %s, Speaker: %s, Date: %s, Token: %s\n',hd{6},hd{7},hd{9},hd{10});
        lst=ffx{3};
        for it=2:size(lst,1);
            [yy,fs,hd,ffx]=readsfs(ffx,lst(it,1),lst(it,2),mode,0);
            nf=hd{1}(7);
            fd=hd{1}(9);
            fprintf(1,'%3d.%02d %ss @ %sHz (%d frames): %s\n',lst(it,1),lst(it,2),sprintsi(nf*fd),sprintsi(1/fd),nf,ffx{4}{it,1});
        end
    end
else
    it=[];
    xfid=[];                  % xfid will be non-empty second time around
    while (isempty(it))                 % may go round this loop twice
        if ischar(ff)           % If ff is a string we must read file
            if any(mode=='d')
                ff=fullfile(voicebox('dir_data'),ff);
            end
            fid=fopen(ff,'rb','b');
            if fid == -1 error('Can''t open %s for input',ff); end
            
            t=fread(fid,512,'uchar').';
            if (t(1:3)~='UC2')
                error(sprintf('%s is not an SFS file type UC2',ff))
            end
            itemlist = [0 1 0 0 t(512)];
            proglist={};
            for i=2:200
                pos = ftell(fid);
                [t,n]=fread(fid,512,'uchar');
                if (n < 512) break; end
                mm=pow2(1,8*([0 1 2 3]+(t(512)==0)*[3 1 -1 -3]));
                itemlist(i,:)=[mm*[t(389:392) t(393:396) t(413:416)] pos t(512)];
                if itemlist(i,1)>=2^31
                    itemlist(i,1)=itemlist(i,1)-2^32;
                end
                if abs(itemlist(i,1))>29
                    if any(mode=='x')
                        itemlist(i,:)=[];
                        break;
                    else
                        error('%s: %d is not a valid SFS item type',ff,itemlist(i,1));
                    end
                end
                proglist{i,1}=char(zerotrim(t(1:256)'));
                proglist{i,2}=char(zerotrim(t(257:384)'));
                proglist{i,3}=char(zerotrim(t(437:456)'));
                seekstat=fseek(fid,itemlist(i,3),'cof');
                if seekstat
                    if any(mode=='x')
                        itemlist(i,:)=[];
                        break;
                    else
                        error('%s: Unexpected end of file',ff);
                    end
                end
            end
            ffx={ff; [fid 0 0 0]; itemlist; proglist};
        else
            ffx=ff;
            ff=ffx{1};
            fid=ffx{2}(1);
            if fid<0
                fid=fopen(ffx{1},'rb',char('b'+(ffx{2}(2)~=0)*('l'-'b')));
            end
        end
        
        % now try to find the requested item
        
        list=ffx{3};
        if nargin<2 ty=0; end
        if nargin<3 sub=0; end
        if ty<0 ty=list(1,1); end
        if ischar(sub)
            lsub=length(sub);
            proglist=ffx{4};
            
            for itt=size(proglist,1):-1:2
                if list(itt,1)==ty & length(proglist{itt,1})>=lsub
                    if strcmpi(sub,proglist{itt,1}(1:lsub))
                        it=itt;
                    end
                end
            end
            if (isempty(it))
                if any(mode=='c') & isempty(xfid)      % try to create item if we haven't tried before
                    xfid=-1;
                    if nargin>=7
                        xname=xpath;
                        xfid=fopen(xname);
                    end
                    if xfid<0
                        if any('/\'==xpath(end))        % would be better to use fullfile()
                            xname=[xpath sub EXESUF];
                        else
                            xname=[xpath '/' sub EXESUF];
                        end
                        xfid=fopen(xname);
                    end
                    if xfid<0
                        error(sprintf('Cannot find executable program %s',sub));
                    else
                        fclose(xfid);
                        fclose(fid); % close this file
                        doscom=['cmd /c "' xname '" ' ffx{1}];
                        %fprintf(1,'Executing: %s\n',doscom);
                        if dos(doscom) % run the program
                            error(sprintf('Error running DOS command: %s',doscom));
                        end
                        ff=ffx{1};          % force reread of header information
                    end
                else
                    error(sprintf('Cannot find item %d.%s in file: %s',ty,sub,ff));
                end
            end
        else % numeric subitem specification
            if sub>0
                it = find(list(:,1)==ty & list(:,2)==sub);
            elseif sub==0
                it = min(find(list(:,1)==ty));
            else
                it = max(find(list(:,1)==ty));
            end
            if (isempty(it))
                error(sprintf('Cannot find item %d.%d in file: %s',ty,sub,ff));
            end
        end
    end % loop up to two times while (isempty(it))
    lit=list(it,:);
    if ffx{2}(3)~=it
        ffx{2}(3)=it;
        ffx{2}(4)=0;
    end
    
    % read the selected item with the correct byte order
    
    if lit(5)~=ffx{2}(2)
        fclose(fid);
        fid=fopen(ffx{1},'rb',char('b'+(lit(5)~=0)*('l'-'b')));
        ffx{2}(1:2)=[fid lit(5)];
        if fid == -1 error(sprintf('Can''t open %s for input',ff)); end
    end
    fseek(fid,lit(4),'bof');
    
    
    y=[];
    fs=0;
    if ~lit(1)                              % read main header
        mb=fread(fid,512,'uchar').';
        if nargout>2
            mc=[1 4; 9 28; 29 32; 37 56; 57 76; 77 96; 97 116; 117 136; 137 296; 297 304; 305 312; 313 332; 333 412];
            hd=cell(14,1);
            hd{1}=[pow2(1,8*([2 3 0 1]+(mb(512)==0)*[1 -1 1 -1]))*[mb(5:8); mb(33:36)].' mb(512)];
            hd{1}(1)=hd{1}(1)/86400+719529;  % convert date format
            for i=1:13
                hd{i+1}=char(zerotrim(mb(mc(i,1):mc(i,2))));
            end
        end
    else
        
        % read the item header
        
        hd=cell(4,1);
        hdr=zeros(1,14);
        
        [str,n]=fread(fid,256,'uchar');
        if (n<256) error(sprintf('Error reading item %d.%d in file: %s',ty,sub,file)); end
        hd{2}=char(zerotrim(str'));
        
        [str,n]=fread(fid,128,'uchar');
        if (n<128) error(sprintf('Error reading item %d.%d in file: %s',ty,sub,file)); end
        hd{3}=char(zerotrim(str'));
        
        hdr(1:8)=fread(fid,8,'long');
        hdr(9)=fread(fid,1,'double');
        if hdr(9) fs=1/hdr(9); end
        hdr(10)=fread(fid,1,'long');
        hdr(11)=fread(fid,1,'double');
        
        [str,n]=fread(fid,20,'uchar');
        if (n<20) error(sprintf('Error reading item %d.%d in file: %s',ty,sub,file)); end
        hd{4}=char(zerotrim(str'));
        
        [hdr(12:14),n]=fread(fid,3,'long');
        if (n<3) error(sprintf('Error reading item %d.%d in file: %s',ty,sub,file)); end
        fseek(fid,44,'cof');
        hd{1}=hdr;
        hd{1}(1)=hd{1}(1)/86400+719529;  % convert date format
        
        % now read the actual data
        
        if nargin<6 nskip=ffx{2}(4);
        elseif nskip<0 nskip=ffx{2}(4);
        end
        
        ksamples=hdr(7)-nskip;
        if nargin>4
            if nmax>=0
                ksamples=min(nmax,ksamples);
            end
        end
        
        if ksamples>0
            ffx{2}(4)=nskip+ksamples;
            fsz=hdr(6);
            if(hdr(10)==1)		% data present
                if(hdr(4)>=0)		% non-structured
                    ds=hdr(5);
                    if(hdr(4)>0)
                        if(ds==4) fmt='float';
                        elseif (ds==8) fmt='double';
                        else error('error in sfs file'); end
                    else
                        if(ds==1 & all(mode~='i')) fmt='uchar';
                        elseif(ds<=2) fmt='short'; fsz=ceil(fsz*ds/2);
                        elseif(ds==4) fmt='long';
                        else error('error in sfs file'); end
                    end
                    fseek(fid,lit(4)+512+nskip*fsz*ds,'bof');
                    nd=fsz*ksamples;
                    [y,n]=fread(fid,nd,fmt);
                    if (n<nd) error(sprintf('Error reading item %d.%d in file: %s',ty,sub,file)); end
                    y = reshape(y,fsz,ksamples)';
                else
                    if (hdr(2)==5)
                        y = cell(ksamples,3);
                        for ifr=1:nskip
                            lf=fread(fid,1,'uchar');
                            fseek(fid,lf,'cof');
                        end
                        
                        for ifr=1:ksamples
                            lf=fread(fid,1,'uchar');
                            tdat=fread(fid,2,'long');
                            y(ifr,:)={tdat(1) tdat(2) char(fread(fid,lf-9,'uchar').')};
                            lf=fread(fid,1,'uchar');
                        end
                    else
                        error(sprintf('Cannot convert item %d.%d in file: %s',ty,sub,file));
                    end
                end
            end
        end
    end
    if all(mode~='f') fclose(fid); ffx{2}(1)=-1; end
end






