function [w,s]=fram2wav(x,tt,mode)
%FRAM2WAV  converts frame values to a continuous waveform [W]=(X,TT,MODE)
%  Inputs:
%          x(nf,p)      is the input signal: one row per frame
%	       tt(nf,3)     specifies the frames. Each row has the form [start_sample end_sample flag]
%                       where flag = 1 for the start of a new spurt.
%                       If tt(:,3) is omitted, a new spurt will be started whenever there is a gap
%                       of more than one between the end of one frame and the beginning or the next.
%                       A new spurt is automatically started if x() = NaN.
%          mode         consists of one or more of the following letters:
%                          z for zero-order hold interpolation (i.e. constant within each frame)
%                          l for linear interpolation within each spurt [default]
%
% Outputs:
%          w(n,p)       contains the interpolated waveforms. Their length is n = tt(nf,2)
%          s(ns,2)      gives the starting and ending sample numbers of each spurt (excluding NaN spurts)
%
%    This routine converts frame-based values to continuous waveforms by performing
%    a chosen method of interpolation. Interpolation is restarted at the beginning of each spurt.

%    Bugs/Suggestions
%      (1)   Additional mode option for cubic interpolation
%      (2)   Additional mode option for interpolation in log domain
%      (3)   Additional mode option for x values being
%            frame averages rather than mid-frame values.

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: fram2wav.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<3
    mode='l';
end
[nf,m]=size(x);
n=round(tt(end,2));
w=repmat(NaN,n,m);
nt=size(tt,2);
ix1=ceil(tt(:,1)); % start of frame sample
ix2=floor(tt(:,2)); % end of frame sample

% determine the start and end of spurts

if nt>2
    ty=tt(:,3)>0;   % frame type set by user
else
    ty=zeros(nf,1);
    ty(2:end)=ix1(2:end)>ix2(1:end-1)+1;    % new spurt whenever a gap
end
ty(1)=1;           % first frame always starts a spurt
ty(isnan(x))=1;    % NaN always ends previous spurt
ty(1+find(isnan(x(1:end-1))))=1; % NaN always forces a new spurt
ty=double(ty);
ty(1:end-1)=ty(1:end-1)+2*ty(2:end);
ty(end)=ty(end)+2;   % last frame always ends a spurtw=repmat(NaN,n,m);  % initialize output to all NaN
nx=ix2-ix1+1;

if any(mode=='z')   % zero-order hold
    for i=1:nf
        if nx(i)
            w(ix1(i):ix2(i),:)=repmat(x(i,:),nx(i),1);
        end
    end
else   % linear interpolation is the default   
    ttm=(tt(:,1)+tt(:,2))/2;    % mid point of frame
    ixm=floor(ttm); % end of first half of frame
    for i=1:nf
        if i==176
            i
        end
        if nx(i)
            tyi=ty(i);
            if tyi==3    % use a zero order hold
                w(ix1(i):ix2(i),:)=repmat(x(i,:),nx(i),1);
            else
                nxm=ixm(i)-ix1(i)+1;
                if nxm
                    if tyi==1    
                        grad=(x(i+1,:)-x(i,:))/(ttm(i+1)-ttm(i));
                    else
                        grad=(x(i,:)-x(i-1,:))/(ttm(i)-ttm(i-1));    
                    end
                    w(ix1(i):ixm(i),:)=repmat(x(i,:),nxm,1)+((ix1(i):ixm(i))'-ttm(i))*grad;
                end
                if nx(i)>nxm
                    if tyi==2
                        grad=(x(i,:)-x(i-1,:))/(ttm(i)-ttm(i-1));
                    else
                        grad=(x(i+1,:)-x(i,:))/(ttm(i+1)-ttm(i));
                    end
                    w(ixm(i)+1:ix2(i),:)=repmat(x(i,:),ix2(i)-ixm(i),1)+((ixm(i)+1:ix2(i))'-ttm(i))*grad;
                end
            end
        end
    end
end

% now sort out the start and end spurt positions

ty(isnan(x))=0;    % Don't count NaN spurts
s=repmat(ix1(bitand(ty,1)>0),1,2);
s(:,2)=ix2(bitand(ty,2)>0);
if ~nargout
    tw=(1:n)';
    for i=size(s,1):-1:2
        j=s(i,1);   % start of new spurt
        tw=[tw(1:j-1); tw(j); tw(j:end)];
        w=[w(1:j-1); NaN; w(j:end)];        % insert a NaN to force a plotting break
    end
    plot(tt(:,1:2)',repmat(x(:)',2,1),'r-+',tw,w,'b-');
end

