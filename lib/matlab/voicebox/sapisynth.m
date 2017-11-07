function [x,fs,txt] = sapisynth(t,m)
%SAPISYNTH  text-to-speech synthesize of text string or matrix [X,FS,TXT]=(T,M)
%
%  Usage:         sapisynth('Hello world');          % Speak text
%                 sapisynth([1 2+3i; -1i 4],'j');    % speak a matrix using 'j' for sqrt(-1)
%          [x,fs]=sapisynth('Hello world','k11');    % save waveform at 11kHz
%                 sapisynth('Hello world','fc');     % use a female child voice if available
%
%  Inputs: t  is either a text string or a matrix
%          m  is a mode string containing one or more of the
%             following options (# denotes an integer):
%
%             'l'   x=a cell array containing a list of talkers
%             'l#'  specify talker # (in the range 1:nvoices)
%             'r#'  speaking rate -10(slow) to +10(fast) [0]
%             'k#'  target sample rate in kHz [22]
%             'o'   audio output [default if no output arguments]
%             'O'   unblocked audio output (may result in simultaneous overlapping sounds)
%             'j'   use 'j' rather than 'i' for complex numbers
%             'm','f' 'c','t','a','s' = Male Female Child, Teen, Adult, Senior
%                       specify any combination in order of priority
%             'v'   autoscale volumne to a peak value of +-1
%             'v#'  set volume (0 to 100) [100]
%             'p#'  set pitch -10 to +10 [0]
%             'n#'  number of digits precision for numeric values [3]
%
% Outputs: x    is the output waveform unless the 'l' option is chosen in
%               which case x is a cell array with one row per available
%               voice containing {'Name' 'Gender' 'Age'} where
%               Gender={Male,Female} and Age={Unknown,Child,Teen,Adult,Senior}
%          fs   is the actual sample frequency
%         txt   gives the actual text sring sent to the synthesiser
%
% The input text string can contain embedded command which are described
% in full at http://msdn.microsoft.com/en-us/library/ms717077(v=vs.85).aspx
% and summarised here:
%
% '... <bookmark mark="xyz"/> ...'               insert a bookmark
% '... <context id="date_mdy"> 03/04/01 </context> ...' specify order of dates
% '... <emph> ... </emph> ...'                   emphasise
% '... <volume level="50"> ... </volume> ...'    change volume level to 50% of full
% '... <partofsp part="noun"> XXX </partofsp> ...'      specify part of speech of XXX: unkown, noun, verb modifier, function, interjection
% '... <pitch absmiddle="-5"> ... </pitch> ...'  change pitch to -5 [0 is default pitch]
% '... <pitch middle="5"> ... </pitch> ...'      add 5 onto the pitch
% '... <pron sym="h eh 1 l ow "/> ...'           insert phoneme string
% '... <rate absspeed="-5"> ... </rate> ...'     change rate to -5 [0 is default rate]
% '... <rate speed="5"> ... </rate> ...'         add 5 onto the rate
% '... <silence msec="500"/> ...'                insert 500 ms of silence
% '... <spell> ... </spell> ...'                 spell out the words
% '... <voice required="Gender=Female;Age!=Child"> ...' specify target voice attributes to be Female non-child
%                                                         Age={Child, Teen, Adult, Senior}, Gender={Male, Female}
%
% Acknowledgement: This function was originally based on tts.m written by Siyi Deng

% Bugs/Suggestions:
%  (1) Allow the speaking of structures and cells
%  (2) Allow a blocking call to sound output and/or a callback procedure and/or a status call
%  (3) Have pitch and/or volume change to emphasise the first entry in a matrix row.
%  (4) extract true frequency from output stream

%      Copyright (C) Mike Brookes 2011
%      Version: $Id: sapisynth.m 4767 2014-06-26 16:35:57Z dmb $
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
persistent vv vvi vvj tsou lsou

% Check that we are on a PC

if ~ispc, error('only works on a PC'); end

% decode the options

if nargin<2
    m='';
end
opts=zeros(52,3); % [exists+number specified, value]
lmode=length(m);
i=1;
while i<=lmode
    if i<lmode  % read a following integer if it exists
        [v,nv,e,ni]=sscanf(m(i+1:end),'%d',1);
    else
        nv=0;
        ni=1;
    end
    k=1+double(lower(m(i)))-'a'+26*(m(i)<'a');
    if k>=1 && k<=52
        opts(k,1)=1+nv;
        if nv
            opts(k,2)=v;
        end
        opts(k,3)=i;  % save position in mode string
    end
    i=i+ni;
end

S=actxserver('SAPI.SpVoice');
V=invoke(S,'GetVoices');  % get a list of voices from the registry
nv=V.Count;
if isempty(vv) || size(vvi,1)~=nv
    vv=cell(nv,3);
    vvi=zeros(nv,6);
    ages={'Senior' 'Adult' 'Teen' 'Child'};
    for i=1:nv
        VI=V.Item(i-1);
        vv{i,1}=VI.GetDescription;
        vv{i,2}=VI.GetAttribute('Gender');
        vvi(i,1)=MatchesAttributes(VI,'Gender=Male');
        vvi(i,2)=MatchesAttributes(VI,'Gender=Female');
        vv{i,3}='Unknown';
        for j=1:length(ages)
            if MatchesAttributes(VI,['Age=' ages{j}])
                vv{i,3}=ages{j};
                vvi(i,2+j)=1;
                break
            end
        end
    end
    vvj=vvi;
    % in the matrix below, the rows and columns are in the order Senior,Adult,Teen,Child.
    % Thus the first row gives the cost of selecting a voice with the wrong age when 'Senior'
    % was requested by the user. A voice of unkown age always scores 0 so entries with negative
    % values are preferred over 'unknown' while those with positive values are not.
    % Diagonal elements of the matrix are ignored (hence set to 0) since correct matches are
    % handled earlier with higher priority.
    vvj(:,3:6)=vvj(:,3:6)*[0 0 1 2; 1 0 2 3; 1 0 0 -1; 1 0 -1 0]'; % fuzzy voice attribute matching
end

% deal with the voice selection options

optv=opts([13 6 19 1 20 3],[3 1 2]);
if opts(12)   % if 'l' option specified - we need to get the voices
    if opts(12)>1
        S.Voice = V.Item(mod(opts(12,2)-1,nv));
    else
        x=vv;
        return
    end
elseif any(optv(:,2))
    optv(:,3)=(1:6)';
    optv=sortrows(optv(optv(:,2)>0,:));  % sort in order of occurrence in mode string
    no=size(optv,1);
    optp=zeros(nv,2*no+1);
    optp(:,end)=(1:nv)'; % lowest priority condition is original rank
    optp(:,1:no)=-vvi(:,optv(:,3));
    optp(:,no+1:2*no)=vvj(:,optv(:,3));
    optp=sortrows(optp);
    S.Voice = V.Item(optp(1,end)-1);
end

% deal with the 'r' option

if opts(18)>1  % 'r' option is specified with a number
    S.Rate=min(max(opts(18,2),-10),10);
end

% deal with the 'v' option

if opts(22)>1  % 'r' option is specified with a number
    S.Volume=min(max(opts(22,2),0),100);
end

% deal with the 'k' option

ff=[11025 12000 16000 22050 24000 32000 44100 48000]; % valid frequencies
if opts(11)>1  % 'k' option is specified with a number
    [v,jf]=min(abs(ff/1000-opts(11,2)));
else
    jf=4;  % default is 16kHz
end
fs=ff(jf);

% deal with the 'n' option

if opts(14)>1  % 'r' option is specified with a number
    prec=opts(14,2);
else
    prec=3;
end

M=actxserver('SAPI.SpMemoryStream');
M.Format.Type = sprintf('SAFT%dkHz16BitMono',fix(fs/1000));
S.AudioOutputStream = M;
if ischar(t)
    txt=t;
else
    txt='';
    if numel(t)
        sgns={' minus ', '', ' plus '};
        sz=size(t);
        w=permute(t,[2 1 3:numel(sz)]);
        sz(1:2)=sz(1)+sz(2)-sz(1:2); % Permute the first two dimensions for reading
        szp=cumprod(sz);
        imch='i'+(opts(10)>0);
        vsep='';
        for i=1:numel(w)
            wr=real(w(i));
            wi=imag(w(i));
            switch((wr~=0)+2*(wi~=0))+4*(abs(wi)==1)
                case {0,1}
                    txt=[txt sprintf('%s%.*g',vsep,prec,wr)];
                case 2
                    txt=[txt sprintf('%s%.*g%c,',vsep,prec,wi,imch)];
                case 3
                    txt=[txt sprintf('%s%.*g%s%.*g%c,',vsep,prec,wr,sgns{2+sign(wi)},prec,abs(wi),imch)];
                case 6
                    if wi>0
                        txt=[txt vsep imch ','];
                    else
                        txt=[txt vsep 'minus ' imch ','];
                    end
                case 7
                    txt=[txt sprintf('%s%.*g%s%c,',vsep,prec,wr,sgns{2+sign(wi)},imch)];
            end
            % could use a <silence msec="???"/> command here
            vsep=[repmat('; ',1,find([0 mod(i,szp)]==0,1,'last')-1) ' '];
        end
    end
end

% deal with the 'p' option

if opts(16)>1  % 'r' option is specified with a number
    txt=[sprintf('<pitch absmiddle="%d"> ',min(max(opts(16,2),-10),10)) txt];
end

invoke(S,'Speak',txt);
x = mod(32768+reshape(double(invoke(M,'GetData')),2,[])'*[1; 256],65536)/32768-1;
delete(M);      % delete output stream
delete(S);      % delete all interfaces

if opts(22)==1 % 'v' option with no argument
    x=x*(1/max(abs(x))); % autoscale
end
if opts(15)>0 || opts(41)>0 || ~nargout % 'o' option for audio output
    while opts(41)==0 && ~isempty(tsou) && toc(tsou)<lsou
    end
    sound(x,fs);
    tsou=tic;   % save time
    lsou=length(x)/fs;
end
