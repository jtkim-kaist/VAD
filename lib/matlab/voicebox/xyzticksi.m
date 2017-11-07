function s=xyzticksi(ax,ah)
%XYZTIXKSI labels an axis of a plot using SI multipliers S=(AX,AH)
%
% This routine is not intended to be called directly. See XTICKSI and YTICKSI.

%	   Copyright (C) Mike Brookes 2009
%      Version: $Id: xyzticksi.m 9302 2017-01-18 16:19:20Z dmb $
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

% Note that "mu" = char(181) assumes Western European encoding
% Bugs:
%   (1) ipan=3 or 4 is not really debugged yet:
%   (2) axis lengths incorrect for 3D graphs
%   (3) should take account of axis shortening due to long labels at the ends
%   (4) should calculate axis orentation from CameraPosition, CameraTarget and CameraUpVector
if nargin<2
    ah=gca;
    if nargin<1
        ax=1;
    end
end
axfield={'XLim' 'YLim' 'ZLim'; 'XTick' 'YTick' 'ZTick'; 'XMinorTick' 'YMinorTick' 'ZMinorTick'; 'XTickLabel' 'YTickLabel' 'ZTickLabel'; 'XScale' 'YScale' 'ZScale'};
tryglobal=nargout>0;
digith=1;    % height of a digit in font units
digitw=0.5;    % width of a digit in font units

prefix={'y','z','a','f','p','n','µ','m','','k','M','G','T','P','E','Z','Y'};
marg=[2 0.5 0.25 0.25];     % gap between labels in font units
ntreq=[3 2 2 1];        % minimum number of labelled ticks required as a function of IPAN
% grid template: each pair is [#steps final-value]. Start=10, end=100
lgridtem={1; [1 20 1 50 1]; [1 20 4]; 9; [2 20 8]; [5 20 2 30 7]; [10 20 5 30 4 50 5]};
ngrid=length(lgridtem);
lgrid=cell(ngrid,1);
agrid=zeros(ngrid,1);       % min and max ratio in decades
% build the actual grid layouts
for i=1:ngrid
    igridtem=[lgridtem{i} 100];
    igrid=zeros(1,sum(igridtem(1:2:end)));
    ntem=length(igridtem)/2;
    k=0;
    tick0=10;
    for j=1:ntem
        nstep=igridtem(2*j-1);
        igrid(k+1:k+nstep)=tick0+(0:nstep-1)*(igridtem(2*j)-tick0)/igridtem(2*j-1);
        k=k+nstep;
        tick0=igridtem(2*j);
    end
    agrid(i)=sum(log10([igrid(2:end) 100]./igrid).^2); % average log interval
    lgrid{i}=igrid';                 % grid positions
end
minsubsp=1;        % minimum subtick spacing (in font units)
delcheck=[log10(2) log10(5) 2];     % check linear spacing
delval=[1 2 5];
dosubtick=0;         % default is not to allow subticks

ngrid=length(lgrid);
loggrid=cell(ngrid,1);
for i=1:ngrid
    loggrid{i}=log10(lgrid{i})-1;
end

getgca=get(ah);
set(ah,'Units','points','FontUnits','points');
getgcac=get(ah);
set(ah,'Units',getgca.Units,'FontUnits',getgca.FontUnits); % return to original values
if ax==1
    widthc=getgcac.Position(3)/getgcac.FontSize;     % x axis length in font units
    axdir=[1 0];        % unit vector in axis direction
else
    widthc=2*getgcac.Position(4)/getgcac.FontSize;     % y axis length in font units
        axdir=[0 1];        % unit vector in axis direction
end
axdir=max(abs(axdir),1e-10);        % avoid infinity problems
a=getgca.(axfield{1,ax})(1);
b=getgca.(axfield{1,ax})(2);

ntick=0;
tickint=[];                 % integer label
tickdp=[];                  % tick decimal point position
ticksi=[];                  % tick SI multiplier (always a multiple of 3)
subtick=[];                 % sub tick positions
if strcmp(getgca.(axfield{5,ax}),'log')   % log axis
    width10=widthc/log10(b/a); % fount units per decade
    ai3=3*ceil(log10(a)/3);             % lowest SI multiplier
    bi3=3*floor(log10(b)/3);            % highest SI multiplier
    if ai3>=-24 && bi3<=24              % do nothing if outside the SI multiplier range
        % first sort out if we can use a global SI multiplier
        if tryglobal && a>=10^(bi3-1)
            gi=bi3;
            s=prefix{9+gi/3};
            globalsi=1;                % global multiplier
        else
            gi=0;
            globalsi=0;                % disable global multiplier
            s='';

        end
        g=10^gi;
        ag=a/g;
        bg=b/g;
        al=log10(ag);
        bl=log10(bg);
        ai=ceil(al);
        bi=floor(bl);
        ai3=3*ceil(ai/3);
        bi3=3*floor(bi/3);
        for ipan=1:4                    % panic level: 1=spacious, 2=cramped, 3=any two labels, 4=any label
            % first try labelling only exact SI multiples
            margin=marg(ipan);
            incsi=3*ceil(min((2*digitw+margin)/axdir(1),(digith+margin)/axdir(2))/(3*width10));   % SI increment
            switch ipan
                case {1,2}
                    ticksi=incsi*ceil(ai/incsi):incsi:incsi*floor(bi/incsi);
                case {3,4}
                    ticksi=ai3:incsi:bi3;
            end
            ntick=length(ticksi);
            tickint=ones(1,ntick);
            tickdp=zeros(1,ntick);
            if width10>0.25
                ticki=ai:bi;
                subtick=10.^(ticki(ticki~=3*fix(ticki/3)));      % put subticks at all powers of 10;
            end
            if incsi==3         % no point in trying anything else if incsi>3
                ci=floor(al);   % start of first decade that includes the start of the axis
                cibi=ci:bi;     % ennumerate the decades that cover the entire axis
                ndec=bi-ci+1;   % number of decades
                if globalsi
                    siq0=zeros(1,ndec); % no SI multipliers within the labels if using a global multiplier
                else
                    siq0=3*floor((cibi)/3); % determine the SI multiplier for each of the decades
                end
                siw0=siq0~=0;                % width of SI multiplier
                dpq0=max(siq0-cibi+1,1);    % number of decimal places
                for jgrid=1:ngrid
                    igrid=jgrid-(ipan<=2)*(2*jgrid-ngrid-1);
                    lgridi=lgrid{igrid};
                    ngridi=length(lgridi);
                    intq=reshape(repmat(lgridi,1,ndec).*repmat(10.^(cibi+dpq0-siq0-1),ngridi,1),1,[]);
                    dpq=reshape(repmat(dpq0,ngridi,1),1,[]);
                    msk=dpq>0 & rem(intq,10)==0;
                    intq(msk)=intq(msk)/10;
                    dpq(msk)=dpq(msk)-1;
                    widq=1+floor(log10(intq));
                    widq=digitw*(widq+(dpq>0).*max(1,dpq+2-widq)+reshape(repmat(siw0,ngridi,1),1,[]));
                    logvq=reshape(repmat(loggrid{igrid},1,ndec)+repmat(ci:ndec+ci-1,ngridi,1),1,[]);
                    % mask out any ticks outside the axis range
                    msk=logvq>=al & logvq<=bl;
                    widq=widq(msk);
                    logvq=logvq(msk);
                    % debug1=[intq(msk); -1 min((widq(1:end-1)+widq(2:end)+2*margin)/axdir(1),2*(digith+margin)/axdir(2))<=2*width10*(logvq(2:end)-logvq(1:end-1))];
                    if numel(widq)>=ntreq(ipan) && all(min((widq(1:end-1)+widq(2:end)+2*margin)/axdir(1),2*(digith+margin)/axdir(2))<=2*width10*(logvq(2:end)-logvq(1:end-1)))
                        % success: we have an acceptable pattern
                        ntick=numel(widq);       % number of ticks
                        tickint=intq(msk);      % integer label value (i.e. neglecting decimal point)
                        tickdp=dpq(msk);        % number of decimal places
                        siq=reshape(repmat(siq0,ngridi,1),1,[]);    % SI mltiplier power
                        ticksi=siq(msk);        % SI multiplier power masked
                        subtick=[];             % recalculate any subticks
                        dosubtick=igrid>1;
                        break;                  % do not try any more grid patterns
                    end
                end % alternative grid patterns
                % finally just try a linear increment
                if ntick<5
                    ldeltamin=log10(bg- bg*10^(-min((digitw+margin)/axdir(1),(digith+margin)/axdir(2))/width10));  % smallest conceivable increment
                    ildelta=floor(ldeltamin);
                    ix=find(ldeltamin-ildelta<=delcheck,1);
                    jx=ildelta*3+ix;
                    while 1
                        deltax=floor(jx/3);
                        deltav=delval(jx-3*deltax+1);
                        delta=deltav*10^deltax;
                        multq=ceil(ag/delta):floor(bg/delta);   % multiples of delta to display
                        ntickq=numel(multq);
                        if ntickq<=ntick || ntickq<ntreq(ipan)   % give up now
                            break;
                        end
                        intq=deltav*multq;
                        lintq=floor(log10(intq));
                        siq=3*floor((lintq+deltax)/3);      % SI multiplier
                        dpq=siq-deltax;
                        msk=dpq<0;
                        intq(msk)=intq(msk).*10.^(-dpq(msk));
                        dpq(msk)=0;
                        msk=rem(intq,10)==0 & dpq>0;
                        while any(msk)      % remove unwanted trailing zeros
                            dpq(msk)=dpq(msk)-1;
                            intq(msk)=intq(msk)/10;
                            msk=rem(intq,10)==0 & dpq>0;
                        end
                        widq=1+floor(log10(intq));
                        widq=digitw*(widq+(dpq>0).*max(1,dpq+2-widq)+(siq~=0));
                        logvq=log10(multq)+log10(deltav)+deltax;
                        % debug2=[intq; widq; -1 (widq(1:end-1)+widq(2:end)+2*margin)<=2*width10*(logvq(2:end)-logvq(1:end-1))];
                        if all(min((widq(1:end-1)+widq(2:end)+2*margin)/axdir(1),2*(digith+margin)/axdir(2))<=2*width10*(logvq(2:end)-logvq(1:end-1)))
                            ntick=ntickq;
                            tickint=intq;
                            tickdp=dpq;
                            ticksi=siq;
                            dosubtick=1;
                            break
                        end
                        jx=jx+1;                            % try next coarser spacing
                    end
                end
            end % try grid patterns since at most one exact SI multiple
            if ntick>=ntreq(ipan)
                break% quit if we have enough labels
            end
        end% try next panic level
    end    % check if within SI range
    if ntick
        tickexp=gi+ticksi-tickdp;
        tickpos=tickint .* 10.^tickexp;
        ratthresh=10^(minsubsp/width10);   % min subtick ratio
        if dosubtick       % check for subticks
            subtick=[];
            if ntick>1      % at least two labelled ticks
                stepexp=min(tickexp(1:end-1),tickexp(2:end))-1;
                stepint=round((tickpos(2:end)-tickpos(1:end-1)).*10.^(-stepexp));  % always a multiple of 10
                stepleft=tickint(1:end-1).*10.^(tickexp(1:end-1)-stepexp); % leftmost label in units of 10^stepexp
                subbase=10.^ceil(log10(stepint)-1); % base sub-tick interval in units of 10^stepexp
                substep=[-1 -3 5]*((1+[1; 2; 5]*(subbase./stepleft))>ratthresh); % actual step is 1,2 or 5 times subbase
                substep(stepint~=10*substep)=max(2-substep(stepint~=10*substep),0); % but only >1 if stepint==10
                substep=substep.*subbase; % subtick step in units of 10^stepexp
                for i=1:ntick-1
                    ss=substep(i);
                    sl=stepleft(i);
                    if ss
                        subtick=[subtick (sl+(ss:ss:stepint(i)-ss))*10^stepexp(i)];
                        if i==1 && sl/(sl-ss)>ratthresh
                            subtick=[subtick (sl-(ss:ss:floor((tickpos(1)-a)/(ss*10^stepexp(i)))*ss))*10^stepexp(i)];
                        elseif i==ntick-1 && (1+ss/(sl+stepint(1)))>ratthresh
                            subtick=[subtick (sl+stepint(i)+(ss:ss:floor((b-tickpos(end))/(ss*10^stepexp(i)))*ss))*10^stepexp(i)];
                        end
                    end
                end
            end
        end % if subtick
        [tps,ix]=sort([tickpos subtick]);
        nticks=length(tps);
        ticklab=cell(nticks,1);
        for j=1:nticks
            i=ix(j);
            if i>ntick
                ticklab{j}='';
            else
                ticklab{j}=sprintf(sprintf('%%.%df%%s',tickdp(i)),tickint(i)*10^(-tickdp(i)),prefix{ticksi(i)/3+9});
            end
        end
        if width10<2.5
            set(ah,axfield{3,ax},'off');
        end
        set(ah,axfield{2,ax},tps);
        set(ah,axfield{4,ax},ticklab);
    end

else                    % linear axis
    for ipan=1:4                    % panic level: 1=spacious, 2=cramped, 3=any two labels, 4=any label
        margin=marg(ipan);
        % select a global multiplier
        if tryglobal
            gi=3*floor(log10(max(abs(a),abs(b)))/3);
            s=prefix{9+gi/3};

        else
            gi=0;
            s='';
        end
        g=10^gi;
        ag=a/g;
        bg=b/g;
        width1=widthc/(bg-ag);                  % width of 1 plot unit in font units
        ldeltamin=log10(min((digitw+margin)/axdir(1),(digith+margin)/axdir(2))/width1);        % log of smallest conceivable increment
        ildelta=floor(ldeltamin);
        ix=find(ldeltamin-ildelta<=delcheck,1);
        jx=ildelta*3+ix;
        while 1 % loop trying increasingly coarse labelling
            deltax=floor(jx/3);         % exponent of label increment
            deltav=delval(jx-3*deltax+1); % mantissa of label increment
            delta=deltav*10^deltax;         % actual label inrement
            multq=ceil(ag/delta):floor(bg/delta);   % multiples of delta to display
            ntickq=numel(multq);
            if ntickq<ntreq(ipan)   % give up now if too few labels
                break;
            end
            intq=deltav*multq;
            lintq=floor(log10(abs(intq)+(intq==0)));
            siq=3*floor((lintq+deltax)/3)*~tryglobal;      % SI multiplier, but only if no global multiplier
            dpq=siq-deltax;
            msk=dpq<0;
            intq(msk)=intq(msk).*10.^(-dpq(msk));
            dpq(msk)=0;
            msk=rem(intq,10)==0 & dpq>0;
            while any(msk)      % remove unwanted trailing zeros
                dpq(msk)=dpq(msk)-1;
                intq(msk)=intq(msk)/10;
                msk=rem(intq,10)==0 & dpq>0;
            end
            widq=1+floor(log10(abs(intq)+(intq==0)));
            widq=digitw*(widq+(dpq>0).*max(1,dpq+2-widq)+(siq~=0).*(intq~=0)+(intq<0)); % calculate width of each label
            % debug2=[intq; widq; digith+margin<=axdir(2)*width1*delta (widq(1:end-1)+widq(2:end)+2*margin)<=2*width1*delta];
            if all((widq(1:end-1)+widq(2:end)+2*margin)<=2*axdir(1)*width1*delta) || (digith+margin<=axdir(2)*width1*delta);
                ntick=ntickq;
                tickint=intq;
                tickdp=dpq;
                ticksi=siq;
                if deltav>1 && width1*delta>0.5*deltav          % put explicit subticks if delta = 2 or 5
                    mults=ceil(ag*deltav/delta):floor(bg*deltav/delta);
                    subtick=(mults(deltav*fix(mults/deltav)~=mults))*delta/deltav;
                else
                    subtick=[];
                end
                break                       % do not try any more coarser spacings
            end
            jx=jx+1;                            % try next coarser spacing
        end
        if ntick>=ntreq(ipan)
            break% quit if we have enough labels
        end
    end
    if ntick
        tickexp=gi+ticksi-tickdp;
        tickpos=tickint .* 10.^tickexp;
        [tps,ix]=sort([tickpos subtick*10^gi]);
        nticks=length(tps);
        ticklab=cell(nticks,1);
        for j=1:nticks
            i=ix(j);
            if i>ntick
                ticklab{j}='';
            else
                ticklab{j}=sprintf(sprintf('%%.%df%%s',tickdp(i)),tickint(i)*10^(-tickdp(i)),prefix{(ticksi(i)/3)*(tickint(i)~=0)+9});
            end
        end
        set(ah,axfield{2,ax},tps);
        set(ah,axfield{4,ax},ticklab);
        set(ah,axfield{3,ax},'on');
    end
end
