function [k,v]=v_findpeaks(y,m,w,x)
%V_FINDPEAKS finds peaks with optional quadratic interpolation [K,V]=(Y,M,W,X)
%
%  Inputs:  Y(N,1)   is the input signal (does not work with UInt datatype)
%           M        is mode:
%                       'f' include the first sample if a downward initial slope
%                       'l' include the last sample if an upward final slope
%                       'm' return only the maximum peak
%                       'q' performs quadratic interpolation
%                       'v' finds valleys instead of peaks
%           W        is the width tolerance; a peak will be eliminated if there is
%                    a higher peak within +-W. Units are samples or X values
%           X(N,1)   x-axis locations of Y values [default: 1:length(Y)]
%
% Outputs:  K(P,1)   are the positions in X of the peaks in Y (fractional if M='q')
%           V(P,1)   are the peak amplitudes: if M='q' the amplitudes will be interpolated
%                    whereas if M~='q' then V=Y(K).

% Outputs are column vectors regardless of whether Y is row or column.
% If there is a plateau rather than a sharp peak, the routine will place the
% peak in the centre of the plateau. When the W input argument is specified,
% the routine will eliminate the lower of any pair of peaks whose separation
% is <=W; if the peaks have exactly the same height, the second one will be eliminated.
% Unless the 'f' or 'l' options are given, all peak locations satisfy 1<K<N.
%
% If no output arguments are specified, the results will be plotted.
%

%	   Copyright (C) Mike Brookes 2005
%      Version: $Id: v_findpeaks.m 6564 2015-08-16 16:56:40Z dmb $
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

if nargin<2 || ~numel(m)
    m=' ';
elseif nargin>3
    x=x(:); % x must be a column vector
end
ny=length(y);
if any(m=='v')
    y=-y(:);        % invert y if searching for valleys
else
    y=y(:);        % force to be a column vector
end
dx=y(2:end)-y(1:end-1);
r=find(dx>0);
f=find(dx<0);
k=[]; % set defaults
v=[];
if ~isempty(r) && ~isempty(f)    % we must have at least one rise and one fall
    dr=r;
    dr(2:end)=r(2:end)-r(1:end-1);
    rc=ones(ny,1);
    rc(r+1)=1-dr;
    rc(1)=0;
    rs=cumsum(rc); % = time since the last rise
    df=f;
    df(2:end)=f(2:end)-f(1:end-1);
    fc=ones(ny,1);
    fc(f+1)=1-df;
    fc(1)=0;
    fs=cumsum(fc); % = time since the last fall
    rp=repmat(-1,ny,1);
    rp([1; r+1])=[dr-1; ny-r(end)-1];
    rq=cumsum(rp);  % = time to the next rise
    fp=repmat(-1,ny,1);
    fp([1; f+1])=[df-1; ny-f(end)-1];
    fq=cumsum(fp); % = time to the next fall
    k=find((rs<fs) & (fq<rq) & (floor((fq-rs)/2)==0));   % the final term centres peaks within a plateau
    v=y(k);
    if any(m=='q')         % do quadratic interpolation
        if nargin>3
            xm=x(k-1)-x(k);
            xp=x(k+1)-x(k);
            ym=y(k-1)-y(k);
            yp=y(k+1)-y(k);
            d=xm.*xp.*(xm-xp);
            b=0.5*(yp.*xm.^2-ym.*xp.^2);
            a=xm.*yp-xp.*ym;
            j=(a>0);            % j=0 on a plateau
            v(j)=y(k(j))+b(j).^2./(a(j).*d(j));
            k(j)=x(k(j))+b(j)./a(j); % x-axis position of peak
            k(~j)=0.5*(x(k(~j)+fq(k(~j)))+x(k(~j)-rs(k(~j))));    % find the middle of the plateau
        else
            b=0.25*(y(k+1)-y(k-1));
            a=y(k)-2*b-y(k-1);
            j=(a>0);            % j=0 on a plateau
            v(j)=y(k(j))+b(j).^2./a(j);
            k(j)=k(j)+b(j)./a(j);
            k(~j)=k(~j)+(fq(k(~j))-rs(k(~j)))/2;    % add 0.5 to k if plateau has an even width
        end
    elseif nargin>3 % convert to the x-axis using linear interpolation
        k=x(k);
    end
end
% add first and last samples if requested
if ny>1
    if any(m=='f') && y(1)>y(2)
        v=[y(1); v];
        if nargin>3
            k=[x(1); k];
        else
            k=[1; k];
        end
    end
    if any(m=='l') && y(ny-1)<y(ny)
        v=[v; y(ny)];
        if nargin>3
            k=[k; x(ny)];
        else
            k=[k; ny];
        end
    end
    
    % now purge nearby peaks - note that the decision about which peaks to
    % delete is not unique
    
    if any(m=='m')
        [v,iv]=max(v);
        k=k(iv);
    elseif nargin>2 && numel(w)==1 && w>0
        j=find(k(2:end)-k(1:end-1)<=w);
        while any(j)
            j=j+(v(j)>=v(j+1));
            k(j)=[];
            v(j)=[];
            j=find(k(2:end)-k(1:end-1)<=w);
        end
    end
elseif ny>0 && (any(m=='f') || any(m=='l'))
    v=y;
    if nargin>3
        k=x;
    else
        k=1;
    end
end
if any(m=='v')
    v=-v;    % invert peaks if searching for valleys
end

if ~nargout
    if any(m=='v')
        y=-y;    % re-invert y if searching for valleys
        ch='v';
    else
        ch='^';
    end
    plot(1:ny,y,'-',k,v,ch);
end
