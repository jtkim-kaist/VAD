function [v,x,t,m,ze]=quadpeak(z)
%PEAK2DQUAD find quadratically-interpolated peak in a N-D array
%
%  Inputs:  Z(m,n,...)   is the input array (ignoring trailing singleton dimensions)
%                        Note: a row vector will have 2 dimensions
%
% Outputs:  V        is the peak value
%           X(:,1)  is the position of the peak (in fractional subscript values)
%           T        is -1, 0, +1 for maximum, saddle point or minimum
%           M        defines the fitted quadratic: z = [x y ... 1]*M*[x y ... 1]'
%           ZE       the estimated version of Z

%	   Copyright (C) Mike Brookes 2008
%      Version: $Id: quadpeak.m 713 2011-10-16 14:45:43Z dmb $
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

persistent wz a
% first calculate the fixed matrix, a (can be stored if sz is constant)
sz=size(z);         % size of input array
psz=prod(sz);       % number of elements in input array
dz=numel(sz);       % number of input dimensions
mz=find(sz>1);      % non-singleton dimension indices
nm=numel(mz);       % number of non-singleton dimensions
vz=sz(mz);          % size of squeezed input array
dx=max(mz);         % number of output dimensions
if ~nm              % if the input array is a scalar
    error('Cannot find peak of a scalar');
end
nc=(nm+1)*(nm+2)/2;  % number of columns in A matrix
if min(vz)<3
    error('Need at least 3 points in each non-singleton dimension');
end
if isempty(wz) || numel(wz)~=numel(vz) || ~all(wz==vz)
    wz=vz;
    a=ones(psz,nc);
    ix=(0:psz-1)';
    for i=1:nm
        jx=floor(ix/sz(mz(i)));
        a(:,i+nc-nm-1)=1+ix-jx*sz(mz(i));
        ix=jx;
        a(:,(i^2-i+2)/2:i*(i+1)/2)=a(:,nc-nm:i+nc-nm-1).*repmat(a(:,i+nc-nm-1),1,i);
    end
    a=(a'*a)\a';        % converts to polynomial coeficients {x^2 xy y^2 x y 1]
end

% now find the peak

c=a*z(:);   % polynomial coefficients for this data
w=zeros(nm+1,nm+1);
i=1:(nm+1)*(nm+2)/2;
j=floor((sqrt(8*i-7)-1)/2);
w(i+j.*(2*nm+1-j)/2)=c;
w=(w+w.')/2; % make it symmetrical
mr=w(1:nm,1:nm);
we=w(1:nm,nm+1);
y=-(mr\we);
v=y'*we+w(nm+1,nm+1);  % value at peak

% insert singleton dimensions into outputs

x=zeros(dx,1);
x(mz)=y;
m=zeros(dx+1,dx+1);
mz(nm+1)=dx+1;
m(mz,mz)=w;
if nargout>2
    ev=eig(mr);
    t=all(ev>0)-all(ev<0);
end
if nargout>4
    ze=zeros(sz);
    scp=cumprod([1 sz(1:end-1)]);
    ivec=fix(repmat((0:psz-1)',1,dz)./repmat(scp,psz,1));
    xe=[1+ivec-repmat(sz,psz,1).*fix(ivec./repmat(sz,psz,1)) ones(psz,1)];
    ze=reshape(sum((xe*m).*xe,2),sz);
end


if ~nargout && nm<=2
    % plot the data
    desc={'Maximum','Saddle Point','Minimum'};
    if nargout<=2
        ev=eig(mr);
        t=all(ev>0)-all(ev<0);
    end
    if nm==1
        xax=linspace(1,psz,100);
        plot(xax,c(1)*xax.^2+c(2)*xax+c(3),'-r',1:psz,z(:),'ob',x,v,'^k');
        set(gca,'xlim',[0.9 psz+0.1]);
        ylabel('z');
        xlabel(sprintf('x%d',mz(1)));
        title(sprintf('\\Delta = %s: z(%.2g) = %.2g',desc{t+2},y(1),v));
    else
        ngr=17;
        xax=repmat(linspace(1,vz(1),ngr)',1,ngr);
        yax=repmat(linspace(1,vz(2),ngr),ngr,1);
        zq=(c(1)*xax+c(2)*yax+c(4)).*xax+(c(3)*yax+c(5)).*yax+c(6);
        hold off
        mesh(xax,yax,zq,'EdgeColor','r');
        hold on
        plot3(repmat((1:vz(1))',1,vz(2)),repmat(1:vz(2),vz(1),1),reshape(z,vz),'ob',y(1),y(2),v,'^k');
        hold off
        set(gca,'xlim',[0.9 vz(1)+0.1],'ylim',[0.9 vz(2)+0.1]);
        xlabel(sprintf('x%d',mz(1)));
        ylabel(sprintf('x%d',mz(2)));
        zlabel('z');
        title(sprintf('\\Delta = %s: z(%.2g,%.2g) = %.2g',desc{t+2},y(1),y(2),v));
    end
end


