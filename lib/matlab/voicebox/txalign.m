function [kx,ky,nxy,mxy,sxy]=txalign(x,y,maxt,nsd)
%TXALIGN Find the best alignment of two sets of time markers [KX,KY,N,M,S]=(X,Y,MAXT)
% x and y vectors contain a list of non-decreasing time values
% we find the alignment between them to minimize (x-y)^2 with a penalty of maxt^2
% for every unmatched pair of entries from x and y.
% an x value can only match to the nearest y value in either the +ve or -ve direction
% If nsd is specified then the threshold is defined circularly to be nsd 
% times the std deviation of the matches away from the mean

% e.g. txalign([1 11 21 27 31 42 51],[2 12 15 22 25 32 41 52 61],1.1);

%
% Outputs: kx is a column vector the same length as x. kx(i)=j if
%             x(i) is matched to y(j). kx(i)=0 if x(i) is unmatched.
%          ky is a column vector the same length as y. ky(j)=i if
%             y(j) is matched to x(i). kx(j)=0 if y(j) is unmatched.
%          nxy is the number of matched pairs
%          mxy is the mean of y(j)-x(i) for the matched pairs
%          sxy is the standard deviation (biassed) of y(j)-x(i) for matched pairs
%
% If no outputs are specified, txalign plots a graph showing y-x against mean(x,y)

%      Copyright (C) Mike Brookes 2002
%      Version: $Id: txalign.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html%
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

% d(*,i): i=1 cumulative cost if not the second of a pair
%           2 cumulative cost if the second of a pair
%           3 +1 if second of pair given that next is not second of pair
%             -1 if not allowed to be second of pair
%           4 0 if this is an x value, 1 if a y value
lx=length(x);
ly=length(y);

if nargin<4 % no nsd specified
x(lx+1)=2*abs(y(ly))+1;
y(ly+1)=2*abs(x(lx))+1;
lxy=lx+ly+1;
d=zeros(lxy,4);

if lx>0 & ly>0
c0=maxt^2;
ix=1;
iy=1;
d(1,:)=[0 0 -1 y(1)<x(1)];
vp=0;
for id=2:lxy
    if y(iy)<x(ix) % do y next
        v=y(iy);
        d(id,4)=1;
        iy=iy+1;
    else % do x next
        v=x(ix);
        ix=ix+1;
    end
    if d(id,4)==d(id-1,4)
        d(id,3)=-1;
    else
        d(id,2)=d(id-1,1)-c0+(v-vp)^2;
    end
    if ~d(id-1,3) & d(id-1,1)>=d(id-1,2)
        d(id,1)=d(id-1,2);
        d(id-1,3)=1;
    else
        d(id,1)=d(id-1,1);
    end
    vp=v;
end
if ~d(lxy,3) & d(lxy,1)>=d(lxy,2)
    d(lxy,3)=1;
end

% now do the traceback

ix=lx;
iy=ly;
nxy=0;
mxy=0;
sxy=0;
kx=zeros(lx,1);
ky=zeros(ly,1);
while (ix>0) & (iy>0)
    id=ix+iy+1;
    if d(id,3)>0
        ky(iy)=ix;  % iy aligned with ix
        kx(ix)=iy;
        dt=y(iy)-x(ix);
        nxy=nxy+1;
        mxy=mxy+dt;
        sxy=sxy+dt^2;
        ix=ix-1;
        iy=iy-1;
    else
        ix = ix-1+d(id,4);
        iy = iy-d(id,4);
    end
end
if nxy
    mxy=mxy/nxy;
    sxy=sqrt(sxy/nxy-mxy^2);
end
% eliminate junk and plot results
if ~nargout
    x(end)=[];
    y(end)=[];
    x=x(:);
    y=y(:);
    p=zeros(lxy-nxy-1,3);
    p=[x x -(kx==0); y(ky==0) y(ky==0) ones(ly-nxy,1)];
    p(kx>0,2)=y(kx(kx>0));
    p=sortrows(p);
    p(:,1:2)=p(:,1:2)*[0.5 -1; 0.5 1];
    plot(p(p(:,3)==0,1),p(p(:,3)==0,2),'-*',p(p(:,3)<0,1),p(p(:,3)<0,2),'x',p(p(:,3)>0,1),p(p(:,3)>0,2),'o');
    xlabel('Mean of x and y (x,o = unmatched x,y)');
    ylabel('y-x difference');
end
else % nsd specifies how many std deviations to exclude
    [kx,ky,nxy,mxy,sxy]=txalign(x,y,maxt);
    nxy0=nxy+1;
    while nxy<nxy0
        nxy0=nxy;
        mxy0=mxy;
        [kx,ky,nxy,mxy,sxy]=txalign(x+mxy0,y,nsd*sxy);
    end
    mxy=mxy+mxy0;
    if ~nargout
        txalign(x,y,nsd*sxy); % plot it
    end
end
else
    kx=zeros(lx,1);
    ky=zeros(ly,1);
    nxy=0;
    mxy=0;
    sxy=0;
end