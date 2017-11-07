function [y,mm]=momfilt(x,r,w,m)
%MOMFILT calculates moments of a signal using a sliding window Y=(X,R,W,M)
%
% Inputs: x    is the input signal
%         r    is a list of moments to calculate
%              (+ve = relative to mean, -ve = relative to zero)
%         w    is the window (or just the length of a Hamming window)
%              Note: If the window is asymmetric, you should be aware that it gets
%              flipped in the convolution process
%         m    is the sample of w to use as the centre [default=ceil(length(w)/2+0.5)]
%
%         mm   the actual value of m used. Output point y(i) is based on x(i+m-w:i+m-1).
%
% Example:
%  To calculate a running kurtosis using a Hamming window of length 30:
%             y=momfilt(x,[2 4],30); k=y(:,2)./y(:,1).^2

%	   Copyright (C) Mike Brookes 2007
%      Version: $Id: momfilt.m 713 2011-10-16 14:45:43Z dmb $
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
if nargin < 3
   w=hamming(length(x));
elseif  length(w)==1
   w=hamming(w);
end
lw=length(w);
w=w(:);                 % force to a column vector
if nargin < 4
   m=(1+lw)/2;
end
m=max(round(m),1);
mm=m;
r=round(r(:))';             % force integer row vector of moments

lx=prod(size(x));
lxx=lx+m-1;
xx=zeros(lxx,1);
xx(1:lx)=x(:);              % extend with zeros so filter() works correctly

cw=cumsum(w);
sw=cw(end);
y0=repmat(sw,lxx,1);
lxw=min(lxx,lw);
y0(1:lxw)=cw(1:lxw);
y0(lx+1:lx+m-1)=y0(lx+1:lx+m-1)-cw(1:m-1);      % equivalent to y0=filter(w,1,xx^0);
yd=y0(m:end);
yd(abs(yd)<eps)=1;

nr=length(r);
wlx=ones(lx,1);
wlxx=ones(lxx,1);
y=zeros(lx,nr);
mr=max(abs(r));                         % max moment to calculate
mk=zeros(1,mr);                       % list of moments
mk(-r(r<0))=1;                         % choose the moments we need to calculate
maxr=max(r);
if maxr>1
    mk(1:maxr)=1;
end
ml=find(mk>0);
lml=length(ml);
if lml
    mlx=mk;
    mlx(ml)=1:lml;                        % mapping from moment into ml
    wlml=ones(1,lml);
    xm=filter(w,1,xx(:,wlml).^ml(wlxx,:));    % calculate all the moments
    xm=xm(m:end,:)./yd(:,wlml);                             % remove the useless start values and normalize
end
fr=find(r<0);
if length(fr)
    y(:,fr)=xm(:,mlx(-r(fr)));    % zero-centred moments
end
fr=find(r==0);                              % 0'th moment
if length(fr)
    y(:,fr)=1;
end
fr=find(r==1);                              % 1'st moment about mean
if length(fr)
    y(:,fr)=0;
end
fr=find(r==2);                              % 2'nd moment about mean
if length(fr)
    yfr=xm(:,2)-xm(:,1).^2;
    y(:,fr)=yfr(:,ones(1,length(fr)));
end
if maxr>2
    mon=[1 -1];
    bc=[1 -2 1];
    am=zeros(lx,maxr);
    am(:,1)=xm(:,1);                            % copy the mean across
    ml=2:maxr;
    wlml=ones(1,length(ml));
    am(:,2:end)=xm(:,wlml).^ml(ones(lx,1),:);              % calculate powers of the mean
    for i=3:maxr
        bc=conv(bc,mon);                        % calculate binomial coefficients
        fr=find(r==i);
        if length(fr)
            yfr=xm(:,i)+sum(xm(:,i-1:-1:2).*am(:,1:i-2).*bc(wlx,2:i-1),2)+am(:,i)*(bc(i)+bc(i+1));
            y(:,fr)=yfr(:,ones(1,length(fr)));
        end
    end
end