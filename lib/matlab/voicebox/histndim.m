function [v,t]=histndim(x,b,mode)
%HISTNDIM - generates and/or plots an n-dimensional histogram
%
%  Inputs:  X(m,d)   is the input data: each row is one d-dimensiona data point
%           B(3,d)   specifies the histogram bins.
%                         B(1,:) gives the number of bins in each dmension [default 10]
%                         B(2,:) gives the minimum of the first bin in each dimension [default min(X)]
%                         B(3,:) gives the maximum of the last bin in each dimension [default max(X)]
%                    If B has only one column, the same values are use for al dimensions 
%                    If B(1,i)=0 then that dimension will be ignored (and excluded from V)
%           MODE     is a character string containing a combination of the following:
%                        'z' for zero base in the 2D plot [default base = min(V)]
%                        'p' to scale V as probabilities [default actual counts]
%                        'h' to plot a histogam even if output arguments are present
%
% Outputs:  V        d-dimensional array containing the histogram counts
%           T        d-element cell array. d{i} contains the bin boundary values for
%                    the i'th dimension. The length of d{i} is one more than the number of bins
%                    in that dimension.
%
%                    Note that if any of B(1,:) are zero then the number of dimensions in V and elements
%                    of T will be correspondingly reduced.
%
% Example: histndim(randn(100000,2),[20 -3 3]','pz');

%	   Copyright (C) Mike Brookes 2004
%      Version: $Id: histndim.m 713 2011-10-16 14:45:43Z dmb $
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

[n,d]=size(x);
if nargin<3
    mode=' ';
    if(nargin<2)
        b=repmat(10,1,d);
    end
end

if size(b,2)==1
    b=repmat(b,1,d);
end
if size(b,1)<3
    mi=min(x,[],1);
    ma=max(x,[],1);
    w=(ma-mi)./(b(1,:)-0.001);  % nudge slightly to make sure al points included
    b(3,:)=ma+0.0005*w;
    b(2,:)=mi-0.0005*w;
end

acd=find(b(1,:)>0);
sv=b(1,acd);
nbt=prod(sv);
t=cell(length(acd),1);

% loop through each dimension
k=1;        % indexing factor
ok=repmat(1>0,n,1);
ix=repmat(nbt-sum(cumprod(sv)),n,1);
for i=1:length(acd)
    j=acd(i);
    bw=b(1,j)/(b(3,j)-b(2,j));
    bi=ceil((x(:,j)-b(2,j))*bw);
    ok=ok & (bi>0) & (bi<=b(1,j));
    ix(ok)=ix(ok)+k*bi(ok);
    k=k*b(1,j);
    t{i}=b(2,j)+(0:b(1,j))/bw;
end
v=full(sparse(ix(ok),1,1,nbt,1));
if length(sv)>1
    v=reshape(v,sv);
end
if any(mode=='p')
    v=v/n;
end

if ~nargout | any(mode=='h')
    svg=find(sv>1);
    if length(svg)==1
        j=acd(svg);
        bar(b(2,j)+(0.5:sv(svg)-0.5)*(b(3,j)-b(2,j))/b(1,j),v(:));
    elseif length(svg)==2
        j=acd(svg(1));
        k=acd(svg(2));
        bj=b(1,j);
        bk=b(1,k);
        %     imagesc(b(2:3,k),b(2:3,j),reshape(v,b(1,j),b(1,k)));
        vda=kron(reshape(v,bj,bk),[1 1;1 1]);
        if any(mode=='z')
            ba=0;
        else
            ba=min(vda(:));
        end
        vda=[repmat(ba,1,2*bk+2);repmat(ba,2*bj,1) vda repmat(ba,2*bj,1);repmat(ba,1,2*bk+2)];
        jda=kron(t{svg(1)},[1 1]);
        jda=jda-(jda(3)-jda(2))*0.01*[-0.5 (-1).^(1:2*bj) 0.5]; % nudge slightly to avoid MATLAB plotting bug
        kda=kron(t{svg(2)},[1 1]);
        kda=kda-(kda(3)-kda(2))*0.01*[-0.5 (-1).^(1:2*bk) 0.5];
        surf(jda,kda,vda');
        ylabel(sprintf('Axis %d',k));
        xlabel(sprintf('Axis %d',j));
        colorbar;
    else
        fprintf(2,'Error in %s: Cannot plot 3-D histogram\n',mfilename);
    end
end
