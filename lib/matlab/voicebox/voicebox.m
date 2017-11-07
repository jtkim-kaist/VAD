function y=voicebox(f,v)
%VOICEBOX  set global parameters for Voicebox functions Y=(FIELD,VAL)
%
%  Inputs:  F   is a field name
%           V   is a new value for the field
%
% Outputs:  Y   is set equal to the structure of parameters if the
%               f and v inputs are both present or both absent. If only
%               input f is specified, then y is set to the value of the
%               corresponding field or null if it doesn't exist.
%
% You can override the defaults set here by setting the environment variable "VOICEBOX"
% to the path of an m-file that contains lines like "% PP.dir_temp='F:\TEMP';"
%
% This routine contains default values for constants that are used by
% other functions in the VOICEBOX toolbox. Values in the first section below,
% entitled "System-dependent directory paths" should be set as follows:
%
%    PP.dir_temp     directory for storing temporary files
%    PP.dir_data     default directory to preappend to speech data file names
%                    when the "d" option is specified in READWAV etc.
%    PP.shorten      location of SHORTEN executable. SHORTEN is a proprietary file compression
%                    algorithm that is used for some SPHERE-format files. READSPH
%                    will try to call an external decoder if it is asked to
%                    read such a compressed file.
%    PP.sfsbin       location of Speech Filing Sysytem binaries. If the "c" option
%                    is given to READSFS, it will try to create a requested item
%                    if it is not present in the SFS file. This parameter tells it
%                    where to find the SFS executables.
%    PP.sfssuffix    suffix for Speech Filing Sysytem binaries. READSFS uses this paremeter
%                    to create the name of an SFS executable (see PP.sfsbin above).
% Other values defined in this routine are the defaults for specific algorithm constants.
% If you want to change these, please refer to the individual routines for a fuller description.

% Bugs/Suggestions
%    (1)  Could allow a * at the end of F to act as a wildcard and return/print a part structure

%      Copyright (C) Mike Brookes 2003
%      Version: $Id: voicebox.m 9312 2017-01-19 13:19:13Z dmb $
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

persistent PP
if isempty(PP)

    % System-dependent directory paths and constants

    PP.dir_temp='F:\TEMP';                      % directory for storing temporary files
    PP.dir_data='E:\dmb\data\speech';           % default directory to preappend to speech data file names
    PP.shorten='C:\bin\shorten.exe';            % location of shorten executable
    PP.flac='C:\bin\flac.exe';                  % location of flac executable
    PP.sfsbin='F:\Program Files\SFS\Program';   % location of Speech Filing Sysytem binaries
    PP.sfssuffix='.exe';                        % suffix for Speech Filing Sysytem binaries
    PP.memsize=50e6;                            % Maximum amount of temporary memory to use (Bytes)

    % DYPSA glottal closure identifier

    PP.dy_cpfrac=0.3;           % presumed closed phase fraction of larynx cycle
    PP.dy_cproj=0.2;            % cost of projected candidate
    PP.dy_cspurt=-0.45;         % cost of a talkspurt
    PP.dy_dopsp=1;              % Use phase slope projection (1) or not (0)?
    PP.dy_ewdly=0.0008;         % window delay for energy cost function term [~ energy peak delay from closure] (sec)
    PP.dy_ewlen=0.003;          % window length for energy cost function term (sec)
    PP.dy_ewtaper=0.001;        % taper length for energy cost function window (sec)
    PP.dy_fwlen=0.00045;        % window length used to smooth group delay (sec)
    PP.dy_fxmax=500;            % max larynx frequency (Hz)
    PP.dy_fxmin=50;             % min larynx frequency (Hz)
    PP.dy_fxminf=60;            % min larynx frequency (Hz) [used for Frobenius norm only]
    PP.dy_gwlen=0.0030;         % group delay evaluation window length (sec)
    PP.dy_lpcdur=0.020;         % lpc analysis frame length (sec)
    PP.dy_lpcn=2;               % lpc additional poles
    PP.dy_lpcnf=0.001;          % lpc poles per Hz (1/Hz)
    PP.dy_lpcstep=0.010;        % lpc analysis step (sec)
    PP.dy_nbest=5;              % Number of NBest paths to keep
    PP.dy_preemph=50;           % pre-emphasis filter frequency (Hz) (to avoid preemphasis, make this very large)
    PP.dy_spitch=0.2;           % scale factor for pitch deviation cost
    PP.dy_wener=0.3;            % DP energy weighting
    PP.dy_wpitch=0.5;           % DP pitch weighting
    PP.dy_wslope=0.1;           % DP group delay slope weighting
    PP.dy_wxcorr=0.8;           % DP cross correlation weighting
    PP.dy_xwlen=0.01;           % cross-correlation length for waveform similarity (sec)
    
    % now see if an environment variable has been set
    
    vbenv=winenvar('VOICEBOX');
    if exist(vbenv,'file');     % update with locally defined values if defined
        run(vbenv)
    end

    % now check some of the key values for validity

    if exist(PP.dir_temp,'dir')~=7        % check that temp directory exists
        PP.dir_temp = winenvar('temp');     % else use windows temp directory
    end

    [fnp,fnn,fne]=fileparts(mfilename('fullpath'));
    if exist(PP.shorten)~=2        % check that shorten executable exists
        PP.shorten=fullfile(fnp,'shorten.exe'); % next try local directory
        if exist(PP.shorten)~=2        % check if it exists in local directory
            PP.shorten='shorten.exe'; % finally assume it is on the search path
        end
    end

    if exist(PP.flac)~=2        % check that flac executable exists
        PP.flac=fullfile(fnp,'flac.exe'); % next try local directory
        if exist(PP.flac)~=2        % check if it exists in local directory
            PP.flac='flac.exe'; % finally assume it is on the search path
        end
    end

end
if nargin==0
    if nargout==0
        % list all fields
        nn=sort(fieldnames(PP));
        cnn=char(nn);
        fprintf('%d Voicebox parameters:\n',length(nn));

        for i=1:length(nn);
            if ischar(PP.(nn{i}))
                fmt='  %s = %s\n';
            else
                fmt='  %s = %g\n';
            end
            fprintf(fmt,cnn(i,:),PP.(nn{i}));
        end
    else
        y=PP;
    end
elseif nargin==1
    if isfield(PP,f)
        y=PP.(f);
    else
        y=[];
    end
else
    if isfield(PP,f)
        PP.(f)=v;
        y=PP;
    else
        error(sprintf('''%s'' is not a valid voicebox field name',f));
    end
end