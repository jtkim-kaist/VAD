function [v,fx,fx1]=v_ppmvu(x,fsx,m)
%V_PPMVU calculate PPM, VU or EBU level of an audio signal [V,FX,FX1]=(X,FSX,M)
%
% Usage: (1) v=v_ppmvu(x,fs,'a')      % calculate PPM level of signal x with sampling freq fs
%        (2) [v,f]=v_ppmvu(x,fs,'aw') % calculate PPM + fast version as well
%        (3) v=v_ppmvu(x,fs,'au')     % calculate VU level in linear units rather than dB
%        (4) v=v_ppmvu(x,fs,'c')      % calculate VU level
%        (5) v=v_ppmvu(x,fs,'e')      % calculate EBU loudness level
%        (6) [v1,fx]=v_ppmvu(x1,fs,'a')   
%                 v2=v_ppmvu(x2,fx)   % process in chunks: same as  v_ppmvu([x1; x2],fs,'a')
%
%  Inputs: x = input signal, one **column** per channel
%          fsx = sample frequency or fx output from a previous call
%          m = either a structure with algorithm parameters (see below)
%              or an attack/decay time constsant or a character string:
%              'a' UK PPM characteristic [default]
%              'b' DIN PPM
%              'c' VU American
%              'd' VU French
%              'e' EBU short   [default toggle: 'q']
%              'f' EBU medium  [default toggle: 'q']
%            followed by any combination of modifier toggles:
%              'z' remove mean (not yet implemented)
%              'p' preemphasis (not yet implemented)
%              's' set rise time to zero
%              'o' oversample x 4
%              'w' give fast output (with 0 decay time) as well as slow output
%              'v' plot graph
%              'q' average squared signal
%              'u' output magnitude (or mean square value if 'q') instead of dB
%
% Outputs:  y = selected output (same size as x)
%           fx = cell array holding algorithm state
%                or, if 'w' option specified, the fast version of the output
%           fx1 = cell array holding algorithm state (only if 'w' specified)
%
% Algorithm Parameters:
%           mm   text string with options
%           ta   attach time constant or zero if no attack smoothing (seconds)
%           tm   averaging filter duration or zero (seconds)
%           td   decay time constant or zero if no decay smoothing (seconds)

%      Copyright (C) Mike Brookes 2013
%      Version: $Id: v_ppmvu.m 3387 2013-08-23 12:32:47Z dmb $
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

% define persistent constants
persistent m4 k4 h1 g1 tamd deft p0
if isempty(m4)
    % upsample by 4 using two cascaded half-band filters
    m4=5;   % one quarter of the first filter order
    h1f=kaiser(4*m4+1,5.3)'.*sinc(-m4:0.5:m4);
    h1=h1f(2:2:4*m4);
    k4=3;
    g1f=kaiser(4*k4+1,5.3)'.*sinc(-k4:0.5:k4);
    g1=g1f(2:2:4*k4);
    tamd=[0.00632 0 1.01335; % 'a' UK PPM characteristic [default]
        0.00316 0 0.73830; % 'b' DIN PPM
        0.10433 0 0; % 'c' VU American
        0.13089 0 0; % 'd' VU French
        0 0.4 0; % 'e' EBU momentary
        0 3 0]; % 'f' EBU short
    deft={'','','','','q','q'};
    p0.mm='a';
    p0.ta=0;
    p0.tm=0;
    p0.td=0;
end
% decode input arguments
if iscell(fsx)
    fx=fsx;
    [nx,nc]=size(x);
    fsp=fx{1};
    fsx=fsp(1);
    fy=fsp(2);
    nt=fsp(3);
    nm=fsp(4);
    p=fx{2};
    mop=fx{3};
else
    p=p0;  % default values
    if nargin<3 || isempty(m)
        m='a';
    end
    if ischar(m)
        p.mm=m;
        ix=m(1)-'a'+1;
        if ix<=0 || ix>6
            ix=1;
        end
        p.ta=tamd(ix,1);
        p.tm=tamd(ix,2);
        p.td=tamd(ix,3);
        p.mm=[p.mm deft{ix}]; % add in default toggles
    elseif isstruct(m)
        p=m;
    else
        p.mm=' ';
        p.ta=m;
        p.tm=0;
        p.td=m;
    end
    [nx,nc]=size(x);
    mop=rem(sum(repmat('zpsowvqu',length(p.mm),1)-p.mm(ones(8,1),:)'==0,1),2); % set option toggles
    fy=fsx*(1+3*mop(4)); % effective sample rate
    nt=round(fy*0.1);       % 0.1 s resolution in samples fro MA filter
    nm=round(fy*p.tm/nt);   % number of blocks for MA filter
    fsp=[fsx fy nt nm];
    fx={fsp p mop zeros(m4,nc),zeros(2*m4-1,nc),zeros(k4,nc),zeros(2*k4-1,nc), zeros(1,nc), 0, zeros(nm+1,nc), zeros(0,nc)};
end

% Stage 1: oversampling we use two cascaded half-band filters to upsample by 4

if mop(4)
    nx2=2*nx;
    v=zeros(nx2,nc);
    v(1:2:2*m4,:)=fx{4};
    v(2*m4+1:2:nx2,:)=x(1:nx-m4,:);
    fx{4}=x(nx-m4+1:nx,:);
    [v(2:2:nx2,:),fx{5}]=filter(h1,1,x,fx{5});   % delayed by 2*m4 samples
    nx4=4*nx;
    y=zeros(nx4,nc);
    y(1:2:2*k4,:)=fx{6};
    y(2*k4+1:2:nx4-1,:)=v(1:nx2-k4,:);
    fx{6}=v(nx2-k4+1:nx2,:);
    [y(2:2:nx4,nc),fx{7}]=filter(g1,1,v,fx{7});   % delayed by 4*m4+2*k4 samples
else
    y=x;
end
ty=1/fy;
% Stage 2 make positive
if mop(7)
    y=y.^2;
else
    y=abs(y);
end
% Stage 2: attack filter
if p.ta>0 && ~mop(3)
    za=exp(-ty/p.ta);
    [y,fx{8}]=filter(1-za,[1 -za],y,fx{8});
end

% Stage 3: moving average filter
if nm>0
    ny=size(y,1);
    y=cumsum(y,1);
    nr=fx{9};
    jj=nt-nr:nt:ny; % end of frame indices
    kk=length(jj);
    if kk>0
        nmkk=nm+kk;
        v=zeros(nmkk,nc);
        v(nm+2:nmkk,:)=y(jj(2:kk),:)-y(jj(1:kk-1),:); % v(:,nc) contains the sum of 0.1 second blocks
        v(nm+1,:)=y(jj(1),:)+fx{10}(2,:);
        fx{10}(2,:)=y(ny,:)-y(jj(kk),:);
        v(2:nm,:)=fx{10}(3:nm+1,:);   % saved values
        fx{10}(3:nm+1,:)=v(2+kk:nmkk,:);
        v=cumsum(v,1);
        v(nm+1:nmkk,:)=v(nm+1:nmkk,:)-v(1:kk,:); % perform MA filter
        v(nm,:)=fx{10}(1,:);                % final MA output from previous chunk
        y=v(nm+floor((nr+1:nr+ny)/nt),:)/(nm*nt);   % copy MA ouptuts into output buffer
        fx{10}(1,:)=v(nmkk,:);                % save the final MA filter output for the next chunk
        fx{9}=ny-jj(kk);                    % update the length of the tail
    else                                    % no completed 100 ms blocks in this chunk
        fx{9}(1)=ny+nr;                     % update the length of the tail
        fx{10}(2,:)=y(ny,:)+fx{10}(2,:);    % update the sum of the tail
        y=fx{10}(ones(ny,1),:);             % set outputs to the final MA output from previous chunk
    end
end
% Stage 4: decay filter
if p.td>0
    [v,vk,fx{11}]=maxfilt(y,exp(-ty/p.td),Inf,1,fx{11});
else
    v=y;
end

% Stage 5: decimate and ouput
if mop(5)    % if 'w' option specified, output y as well
    if mop(4) % if oversampling, we need to decimate
        y=y(1:4:end,:);
        v=v(1:4:end,:);
    end
    if ~mop(8)
        y=(20-10*mop(7))*log10(y);
        v=(20-10*mop(7))*log10(v);
    end
    fx1=fx;
    fx=y;
else
    if mop(4) % if oversampling, we need to decimate
        v=v(1:4:end,:);
    end
    if ~mop(8)
        v=(20-10*mop(7))*log10(v);
    end
end
if mop(6) || nargout==0
    t=(1:length(x))/fsx;
    ax(1)=subplot(211);
    plot(t,x,'-b');
    if mop(7)>0 && mop(8)>0
        v=sqrt(v);
        if mop(5)
            y=sqrt(y);
        end
    end
    ax(2)=subplot(212);
    if mop(5)    % if 'w' option specified, output y as well
        plot(t,v,'-b',t,y,'-r');
    else
        plot(t,v,'-b');
    end
    if ~mop(8)
        set(gca,'ylim',[max(min(v)-1,max(v)-40) max(v)+1]);
        ylabel('dB');
    elseif mop(7)
        ylabel('x_{RMS}');
    else
        ylabel('|x|');
    end
    linkaxes(ax,'x');
    xlabel(['Time (' xticksi 's)']);
end
