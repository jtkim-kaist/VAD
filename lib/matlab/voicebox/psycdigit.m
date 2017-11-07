function [m,v]=psycdigit(proc,r,mode,p,q,xp,noise,fn,dfile,ofile)
%PSYCDIGIT measures psychometric function using TIDIGITS stimuli
%
% Usage:
%         (1)[m,v]=psycdigit([],[],'GMPWrn10',[],[],[],[],[],dfile);
%                       % measure SRT using addditive white noise, repetitions allowed, data in .WAV format
%         (2)[m,v]=psycdigit(@specsub,[],'GMPn10',[],[],[],[],[],dfile);
%                       % compare spectral subtraction with unprocessed noisy speech
%
% Inputs:
%         proc    processing function handle (e.g. @specsub) called as follows:
%                    y=proc(x,fs,i,r)  % process noisy waveform x through model i with parameters r
%                    y=proc(x,fs,snr,i,r) % process clean speech degraded to snr through model i with parameters r
%                    y=proc()   % return comment text string describing algorithm (if 'c' option specified)
%         r       parameters for proc (omitted from call if r=[])
%         mode    string containing options (see below for list)
%         p,q,xp  parameters passed on to psycest or psycestu
%         noise   noise waveform or wav file containing noise
%         fn      noise waveform sample frequency [16000]
%         dfile   path to TIdigits folder
%         ofile   output text file (see below for output file format)
%
% Outputs:
%          m(2,n,3) estimated srt and slope of all models
%                   m(i,n,k): i={1=srt, 2=slope}, n=model number, k={1=mean, 2=mode (MAP), 3=marginal mode}
%          v(3,n)   estimated covariance matrix entries:
%                   [var(srt) cov(srt,slope) var(slope)]' of n'th model
%
% List of options for "mode" input:
%         a       proc adds its own noise
%         b*      base figure number for plotting results [100]
%         c       y=proc() returns comment string
%         e/E*    plot evolving psychometric functions *=1,2,3 for mean, mode, marginal mode (F=after each trial)
%         f/F*    plot psychometric functions *=1,2,3 for mean, mode, marginal mode (F=after each trial)
%         g       prompt with number of digits
%         G       prompt with SNR
%         l*      min token length (in digits)
%         L*      max token length (if different from min)
%         m*      use * external models [default 1]
%         M       include an extra model with no processing
%         n*      *=average number of trials per model
%         p/P     plot pdf (P=after each trial)
%         r       allow repetitions
%         s       respond s to save the noisy stimulus as a .wav file
%         t/T*    taper noise */10 seconds at start and end
%         v/V*    plot srt/slope convergence (V= after each trial)
%                 *=1,2,3 for mean, mode, marginal mode
%         x*      add */10 seconds of noise to the front of the speech before processing
%         X*      truncate */10 seconds from the start of the processed speech
%         W       data is in microsoft WAV format
%         z       omit tokens containing "oh"
%         Z       omit tokens containing "zero"
%
% Output file format
%   Each line starts 'x ' where x is one of the following
%       %  Comment
%       V  File version type
%       O  mode options
%       P  details about proc
%       C  Comment returned by proc
%       M  measurement

% Future mode options:
%        [ d     score as single digits ]
%        [ i/I   plot SRT improvement (I=after each trial) ]
%        [ j*    scaling: 0=autoscale each token, 1=constant speech,2=const noise, 3=const noisy ]
%        [N*     ignore the first * trials ]
%        [o/O    write to output file, O write result of each probe]
%        [ S     save all stimuli as wav files ]
%        [ u     do not normalize the noise level ]
%        [ y*    type of noise ]
%
% Bugs/Suggestions
% (1) Add sounds to indicate error and end of test
% (2) Routine could return a label to replace "SNR" in necessary
% (3) Add an input option to "discard" this sample
% (4) output file should include mode argument + date + IP address + computer name + r arguments + p/q values
% (5) Quote filenames in output file or else never have anything else on the line
% (6) silence detection doesn't work

%      Copyright (C) Mike Brookes 2010
%      Version: $Id: psycdigit.m 8166 2016-07-08 14:46:15Z dmb $
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

%     25+4     m*      use * external models [default 1]
%     26    M       include an extra model with no processing
%      1   a       proc adds its own noise
%      27+5   n*      *=average number of trials per model
%     23+2    l*      min token length (in digits)
%      24+3   L*      max token length (if different from min)
%      7   d       score as single digits
%      37   s       speech-shaped noise instead of white
%       19+1  j*      scaling: 0=autoscale each token, 1=constant speech, 2=const noise, 3=const noisy
%      35   r       allow repetitions
%     41    u       unimodal psychometric function
%     32    P       plot pdf
%     31   p       plot srt convergence
%     13    g       prompt with number of digits
%      14   G       prompt with SNR
%      47+6   x*      add */10 seconds of noise to the front of the speech
%      48+7   X*      truncate */10 from the start of the processed speech
% input parameters
persistent tok mtok tigflgp digitsp
if nargin<10
    ofile='';
    if nargin<9
        dfile='F:\home\dmb\data\old_speech\tidigits';
        %         dfile='Y:\Speech\TIDigits\Disc 1\TIDIGITS';
        if nargin<8
            fn=16000;
            if nargin<7
                noise=[];
                if nargin<6
                    xp=[];
                    if nargin<5
                        q=[];
                        if nargin<4
                            p=[];
                            if nargin<3
                                mode='';
                                if nargin<2
                                    r=[];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% parse the mode string

i=1;
pv=[0 3 3 1 25  0  0 1 1 1 1 100]; % pv default values if option is unspecified
px=[0 3 3 1 25 20 18 1 1 1 1 100]; % pv values if option is given without a value
pvs='jlLmnxXvVfFb';
pf=zeros(1,52);
mval=0;
lmode=length(mode);
while i<=lmode
    if i<lmode  % read a following integer if it exists
        [v,nv,e,ni]=sscanf(mode(i+1:end),'%d',1);
    else
        nv=0;
        ni=1;
    end
    j=find(mode(i)==pvs,1);
    k=1+2*(double(lower(mode(i)))-'a')+(mode(i)<'a');
    if k>=1 && k<=52
        pf(k)=1-pf(k);
        if ~isempty(j)
            if nv==0
                pv(j)=px(j);
            else
                pv(j)=v;
                mval=mval || j==4;  % m has a value specified
            end
        end
    end
    i=i+ni;
end
if isempty(proc)
    pv(4)=0; % not allowed any processed models if not process is specified
end
% derived input parameters

varthr=20;   % variance threshold for "silence"
nxevo=30; % size of evloving pdf

if pf(23)>pf(24) || pv(2)>pv(3)
    pv(3)=pv(2);                % Make max token length reasonable
end
if pf(24)>pf(23)
    pv(2)=1;                    % min =1 if only max is specified
end
zmodel=pv(4)+(pf(26)>0); % total number of models (including null model)
ntrial=zmodel*pv(5);  % total number of trials
if pf(12)>pf(11)    % 'F' specified but not 'f'
    pv(10)=pv(11);  % copy across the average type selector
end
if pf(11)+pf(12)>0 && (pv(10)<1 || pv(10)>3)
    error('Invalid average type for option f/F');
end

% sort out p argument to psycest or psycestu
unimode=pf(41);
if isempty(p)
    p=0.75;     % default recognition rate at threshold
end
if size(p,1)<3-unimode
    p(2-unimode,1)=0.04;     % miss probability
    p(3-unimode,1)=0;     % dummy entry for guess probability
end
if size(p,2)<(zmodel)
    p=p(:,1+mod(0:zmodel-1,size(p,2)));
end

% first sort out the digit samples

tigflg=[pv(2:3) pf(51:52)];
if isempty(tok) || any(tigflg~=tigflgp) || ~strcmp(dfile,digitsp)
    disp('Scanning TIDIGITS files ...');
    digitsp=dfile;
    tigflgp=tigflg;
    if any(dfile(end)=='/\')
        dfile(end)=[];  % remove a trailing separator
    end
    dirlist{1}='';
    ntok=0;
    mtok=zeros(1,pv(3));
    tok=cell(1,2);
    while ~isempty(dirlist)
        dd=dir([dfile dirlist{1}]);
        for i=1:length(dd)
            name=dd(i).name;
            if name(1)~='.'   % ignore directories starting with '.'
                if dd(i).isdir
                    dirlist{end+1}=[dirlist{1} '\' name];
                elseif length(name)>4 && strcmpi(name(end-3:end),'.wav')
                    digs=name(1:end-4);
                    digz=upper(digs)=='Z';
                    digo=upper(digs)=='O';
                    digs(digo | digz)='0';
                    digs=digs(digs>='0' & digs<='9');
                    ndigs=length(digs);
                    if ndigs>=tigflg(1) && ndigs<=tigflg(2) && ~(tigflg(3) && any(digo)) && ~(tigflg(4) && any(digz))
                        ntok=ntok+1;
                        mtok(ndigs)=mtok(ndigs)+1;
                        tok{ntok,1}=[dirlist{1} '\' name];
                        tok{ntok,2}=digs;
                    end
                end
            end
        end
        dirlist(1)=[];  % remove this directory from the list
    end
end
ntok=size(tok,1);
if pf(46)
    [s,fs]=readwav([dfile tok{1,1}]); % get the first speech token to set fs
else
    [s,fs]=readsph([dfile tok{1,1}]); % get the first speech token to set fs
end

% calculate guess probability assuming you get the correct number of digits

if any(mode=='d')
    p(3-unimode,:)=0.1;
else
    p(3-unimode,:)=0.1.^(1:pv(3))*mtok'/ntok;
end

% now initialize the models

if unimode
    [xx,ii,m,v]=psycestu(-zmodel,p,q,xp); % initialize all models
    [pact,qact]=psycestu(0);  % save the actual parameters
else
    [xx,ii,m,v]=psycest(-zmodel,p,q,xp); % initialize all models
    [pact,qact]=psycest(0);  % save the actual parameters
end
if pv(4) && pf(5)
    x=[];  % set x=[] to force a description output
    ii=1;   % for now, only do model 1
    if pf(1)   % if process adds its own noise
        if isempty(r)  % proc does not require any auxilliary parameters
            if mval
                procdesc=proc(x,fs,xx,ii);  % process the noisy speech
            else
                procdesc=proc(x,fs,xx);  % process the noisy speech
            end
        else
            if mval
                procdesc=proc(x,fs,xx,ii,r);  % process the noisy speech
            else
                procdesc=proc(x,fs,xx,r);  % process the noisy speech
            end
        end
    else   % if process does not add its own noise
        if isempty(r)  % proc does not require any auxilliary parameters
            if mval
                procdesc=proc(x,fs,ii);  % process the noisy speech
            else
                procdesc=proc(x,fs);  % process the noisy speech
            end
        else
            if mval
                procdesc=proc(x,fs,ii,r);  % process the noisy speech
            else
                procdesc=proc(x,fs,r);  % process the noisy speech
            end
        end
    end
else
    procdesc='';
end

% now initialize the output file

if pf(29) || pf(30)   % o/O output info
    nw=fix(datevec(now));
    of=sprintf('psy%4d%02d%02d%02d%02d.txt',nw(1:5)); % filename includes date and time to nearest minute
    ofid=fopen(of,'wt');
    if ~ofid
        error('Cannot write to %s',of);
    end
    fprintf(ofid,'%% %s evaluation on %s\n',mfilename,datestr(now));
    fmfnm=[mfilename('fullpath') '.m'];
    dd=dir(fmfnm);
    fprintf(ofid,'%% %s = %d bytes = %s\n',fmfnm,dd.bytes,dd.date);
    fprintf(ofid,'V %d\n',2); % print file format version number
    fmfnm=which(func2str(proc));
    dd=dir(fmfnm);
    fprintf(ofid,'P %s = %d bytes = %s\n',fmfnm,dd.bytes,dd.date);
    if numel(procdesc)
        fprintf(ofid,'C %s\n',procdesc);
    end
    fprintf(ofid,'O %s\n',mode);
end

% now start testing

disp(['Testing ' procdesc]);
% now do the main loop
mnt=zeros(2+2*unimode,zmodel,3,ntrial+1);
vnt=zeros(3+unimode,zmodel,ntrial+1);
mnt(:,:,:,1)=m;
vnt(:,:,1)=v;
i=0;
imax=0;
quitit=0;
while ~quitit
    i=i+1;
    isp=min(1+floor(rand(1)*ntok),ntok); % select a token
    if pf(46)
        [s,fs]=readwav([dfile tok{isp,1}]); % get the speech token
    else
        [s,fs]=readsph([dfile tok{isp,1}]); % get the speech token
    end
    s=[zeros(pv(6)*round(fs/10),1); activlev(s(:),fs,'n')]; % preappend zeros and normalize speech level
    if pf(1) && ii<=pv(4)  % if process adds its own noise
        if isempty(r)
            if mval
                y=proc(s,fs,xx,ii);  % process the noisy speech
            else
                y=proc(s,fs,xx);  % process the noisy speech
            end
        else
            if mval
                y=proc(s,fs,xx,ii,r);  % process the noisy speech
            else
                y=proc(s,fs,xx,r);  % process the noisy speech
            end
        end
    else
        nn=randn(length(s),1);
        x=nn+s*10^(xx/20);   % create the data
        if ii<=pv(4)
            if isempty(r)
                if mval
                    y=proc(x,fs,ii);  % process the noisy speech
                else
                    y=proc(x,fs);  % process the noisy speech
                end
            else
                if mval
                    y=proc(x,fs,ii,r);  % process the noisy speech
                else
                    y=proc(x,fs,r);  % process the noisy speech
                end
            end
        else
            y=x;            % unprocessed for last model
        end
    end
    y=y(1+pv(7)*round(fs/10):end);   % remove junk from the start
    if pf(13)
        prg=sprintf(' %d',length(tok{isp,2}));
    else
        prg='';
    end
    if pf(14)
        prG=sprintf('SNR=%d dB, ',round(xx));
    else
        prG='';
    end
    if pf(35)
        prr=sprintf(', ''r'' to repeat');
    else
        prr='';
    end
    prompt=[prG 'enter' prg ' digits (''q'' to quit' prr '): '];
    %     ansr=-(var(y)>varthr);  % meant to detect silences but doesn't work
    ansr=-1;
    say=1;
    while ansr<0
        if say
            tdel=0;
            tic;
            soundsc(y,fs);      % *** probably shouldn't be ...sc
            say=0;
        end
        rv=input(prompt,'s');
        tdel=toc;
        if ~isempty(rv)
            if lower(rv(1))=='q'
                quitit=1;
                ansr=2;
            elseif lower(rv(1))=='r' && pf(35)
                say=1;
            elseif lower(rv(1))=='s' && pf(37)   % save the token
                ofn=input('File name: ','s');
                if numel(ofn)
                    writewav(y,fs,ofn);
                end
            elseif all(rv>='0') && all(rv<='9') && ( ~pf(13) || length(rv)==length(tok{isp,2}))
                ansr=strcmp(rv,tok{isp,2});
            end
        end
    end
    quitit=quitit || i==ntrial;   % mark as quiting if we have done all the trials
    jj=ii;  % remember which model has just been updated
    xxold=xx; % and what the SNR was
    if ansr<2  % valid answer: update the pdfs
        if unimode
            [xx,ii,m,v]=psycestu(ii,xx,ansr);
        else
            [xx,ii,m,v]=psycest(ii,xx,ansr);
        end
        mnt(:,:,:,i+1)=m;
        vnt(:,:,i+1)=v;
        imax=i;
    end
    if pf(30) || quitit && pf(29)       % 'O/o': output
        if ansr>1
            rv=num2str(double(rv(1)));
        end
        % could add in token name and length
        fprintf(ofid,'M %d %d %.3g %s %s %d %.1f',i,ii,xxold,tok{isp,2},rv,ansr,tdel);
        fprintf(ofid,' %.3g',m(:,jj,:),v(:,jj));
        fprintf(ofid,'\n');
    end
    if pf(32)                           % 'P': plot PDF: figures 1:m
        figure(jj);
        if unimode
            psycestu(jj);
        else
            psycest(jj);
        end
    elseif quitit && pf(31)             % 'p': plot PDF: figures 1:m
        for jj=1:zmodel
            figure(jj);
            if unimode
                psycestu(jj);
            else
                psycest(jj);
            end
        end
    end
    if pf(12) || quitit && pf(11)       % 'F/f': plot Psychometric function on figure 101
        figure(pv(12)+1);
        if unimode
            qqu.pk=m(1,jj,pv(10));      % peak position
            qqu.w=m(2,1,pv(10));        % peak width
            qqu.ll=m(3,1,pv(10));       % peak slope on low side
            qqu.lh=m(4,1,pv(10));       % peak slope on high side
            qqu.gu=pact(2);             % guess rate
            qqu.la=pact(1);             % lapse rate
            sw=qqu.w*qqu.ll*qqu.lh/(qqu.ll+qqu.lh);   % normalized distance of inflections from peak
            xax=linspace(qqu.pk-(4+sw)/qqu.ll,qqu.pk+(4+sw)/qqu.lh,200);
            bs=(qqu.lh-qqu.ll)/2;
            bu=(qqu.lh+qqu.ll)/2;
            plot(xax,qqu.gu+(1-qqu.gu-qqu.la)*(1+2*exp(-sw)*cosh(bs*(xax-qqu.pk)+bu*abs(xax-qqu.pk))+exp(-2*sw)).^(-1));
        else
            sd=(pact(1,:)-pact(3,:)).*(1-pact(2,:)-pact(1,:))./(m(2,:,pv(10)).*(1-pact(3,:)-pact(2,:)));
            md=m(1,:,pv(10))-sd.*log((pact(1,:)-pact(3,:))./(1-pact(2,:)-pact(1,:)));
            xax=linspace(min(md-3*sd),max(md+3*sd),100);
            for jj=1:zmodel
                plot(xax,psychofunc('',[pact(1,jj); m(:,jj,pv(10)); pact(2:3,jj); qact(10)],xax));
                hold on
            end
            hold off
        end
        axis([xax(1) xax(end) 0 1]);
        xlabel('SNR (dB)');
        ylabel('Recognition Probability');
        title(sprintf('Intelligibility: %s',procdesc));
    end
    if pf(10) || quitit && pf(9)        % 'E/e': plot evolving Psychometric function on figure 103
        figure(pv(12)+3);
        psyevo=zeros(ntrial+1,nxevo);   % space for evolving pdf
        if unimode                      % unimodal psychometric function
            qqu.pk=m(1,jj,pv(10));      % peak position
            qqu.w=m(2,1,pv(10));        % peak width
            qqu.ll=m(3,1,pv(10));       % peak slope on low side
            qqu.lh=m(4,1,pv(10));       % peak slope on high side
            qqu.gu=pact(2);             % guess rate
            qqu.la=pact(1);             % lapse rate
            sw=qqu.w*qqu.ll*qqu.lh/(qqu.ll+qqu.lh);   % normalized distance of inflections from peak
            xax=linspace(qqu.pk-(4+sw)/qqu.ll,qqu.pk+(4+sw)/qqu.lh,nxevo);
            for iet=1:imax+1
                qqu.pk=mnt(1,1,pv(10),iet);     % peak position
                qqu.w=mnt(2,1,pv(10),iet);    	% peak width
                qqu.ll=mnt(3,1,pv(10),iet);     % peak slope on low side
                qqu.lh=mnt(4,1,pv(10),iet);     % peak slope on high side
                sw=qqu.w*qqu.ll*qqu.lh/(qqu.ll+qqu.lh);   % normalized distance of inflections from peak
                bs=(qqu.lh-qqu.ll)/2;
                bu=(qqu.lh+qqu.ll)/2;
                psyevo(iet,:)=qqu.gu+(1-qqu.gu-qqu.la)*(1+2*exp(-sw)*cosh(bs*(xax-qqu.pk)+bu*abs(xax-qqu.pk))+exp(-2*sw)).^(-1);
            end
            imagesc(xax,0:ntrial,psyevo);
            axis 'ij'  % put trial 0 at the top
        else                            % monotonic psychometric function (not implemented)
        end
        xlabel('SNR (dB)');
        ylabel('After trial');
        title(sprintf('Intelligibility: %s',procdesc));
    end
    if pf(44) || quitit && pf(43)   % 'V/v': plot convergence on figure 102
        figure(pv(12)+2);
        if unimode                      % unimodal psychometric function
            sw=mnt(2,1,pv(10),1:imax+1).*mnt(3,1,pv(10),1:imax+1).*mnt(4,1,pv(10),1:imax+1)./(mnt(3,1,pv(10),1:imax+1).*mnt(4,1,pv(10),1:imax+1));
            subplot(221)
            plot(0:imax,squeeze(mnt(1,1,pv(10),1:imax+1)));
            set(gca,'xlim',[0 ntrial]);
            xlabel('After trial');
            ylabel('Peak position (dB SNR)');
            subplot(222)
            plot(0:imax,squeeze(mnt(1,1,pv(10),1:imax+1)-sw.*mnt(3,1,pv(10),1:imax+1)),'-b',0:imax,squeeze(mnt(1,1,pv(10),1:imax+1)+sw.*mnt(4,1,pv(10),1:imax+1)),'-b');
            set(gca,'xlim',[0 ntrial]);
            xlabel('After trial');
            ylabel('Inflections (dB)');
            subplot(223)
            plot(0:imax,squeeze(mnt(3,1,pv(10),1:imax+1)));
            set(gca,'xlim',[0 ntrial]);
            xlabel('After trial');
            ylabel('Upwards slope (prob/dB)');
            subplot(224)
            plot(0:imax,squeeze(mnt(4,1,pv(10),1:imax+1)));
            set(gca,'xlim',[0 ntrial]);
            xlabel('After trial');
            ylabel('Downwards slope (prob/dB)');
        else                            % monotonic psychometric function
            subplot(211);
            for jj=1:zmodel
                plot(0:imax,squeeze(mnt(1,jj,pv(10),1:imax+1)));
                hold on
            end
            hold off
            set(gca,'xlim',[0 ntrial]);
            xlabel('After trial');
            ylabel('SRT (dB SNR)');
            subplot(212);
            for jj=1:zmodel
                plot(0:i,squeeze(mnt(2,jj,pv(10),1:imax+1)));
                hold on
            end
            hold off
            set(gca,'xlim',[0 ntrial]);
            xlabel('After trial');
            ylabel('Slope (Prob/dB)');
        end
    end
end   % main loop for each probe value
if pf(29) || pf(30)
    fclose(ofid);
end
