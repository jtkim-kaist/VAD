function [x,g,xn,gg] = kmeanhar(d,k,l,e,x0)
%KMEANS Vector quantisation using K-harmonic means algorithm [X,G,XN,GG]=(D,K,L,E,X0)
%
%  Inputs:
%
%    D(N,P)  contains N data vectors of dimension P
%    K       is number of centres required
%    L       integer portion is max loop count, fractional portion
%            gives stopping threshold as fractional reduction in performance criterion
%    E       is exponent in the cost function. Significantly faster if this is an even integer. [default 4]
%    X0(K,P) are the initial centres (optional)
%            Alternatively, X0 can be a character determining the initialization method:
%                'f'    Initialize with K randomly selected data points [default]
%                'p'    Initialize with centroids and variances of random partitions
%
%  Outputs:
%
%    X(K,P)  is output row vectors
%    G       is the final performance criterion value (normalized by N)
%    XN      nearest centre for each input point
%    GG(L+1) value of performance criterion before each iteration and at end
%
% The k-harmonic means algorithm selects K cluster centres to minimize 
%                           sum_n(K/sum_k((d_n-x_k)^-e))
% where sum_n is over the N inputs points d_n and sum_k is over the K cluster centres x_k.
%
% It is often a good idea to scale the input data so that it has equal variance in each
% dimension before calling KMEANHAR so that approximately equal weight is given
% to each dimension in the distance calculation.

%  [1] Bin Zhang, "Generalized K-Harmonic Means - Boosting in Unsupervised Learning",
%      Hewlett-Packartd Labs, Technical Report HPL-2000-137, 2000 [Zhang2000]
%      http://www.hpl.hp.com/techreports/2000/HPL-2000-137.pdf

%  Bugs:
%      (1) Could use nested blocking to allow very large data arrays
%      (2) Could then allow incremental calling with partial data arrays (but messy)

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: kmeanhar.m 713 2011-10-16 14:45:43Z dmb $
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

% sort out the input arguments

if nargin<5
    x0='f';
    if nargin<4
        e=[];
        if nargin<3
            l=[];
        end
    end
end
if isempty(e)
    e=4;  % default value
end
if isempty(l)
    l=50+1e-3; % default value
end
sd=5;       % number of times we must be below threshold


% split into chunks if there are lots of data points

memsize=voicebox('memsize');
[n,p] = size(d);
nb=min(n,max(1,floor(memsize/(8*p*k))));    % block size for testing data points
nl=ceil(n/nb);                  % number of blocks

% initialize if X0 argument is not supplied

if ischar(x0)
    if k<n
        if any(x0=='p')                  % Initialize using a random partition
            ix=ceil(rand(1,n)*k);       % allocate to random clusters
            ix(rnsubset(k,n))=1:k;      % but force at least one point per cluster
            x=zeros(k,p);
            for i=1:k
                x(i,:)=mean(d(ix==i,:),1);
            end
        else                                % Forgy initialization: choose k random points [default]
            x=d(rnsubset(k,n),:);         % sample k centres without replacement
        end
    else
        x=d(mod((1:k)-1,n)+1,:);    % just include all points several times
    end
else
    x=x0;
end
eh=e/2;
th=l-floor(l);
l=floor(l)+(nargout>1);   % extra loop needed to calculate final performance value
if l<=0
    l=100;      % max number of iterations ever
end
if th==0
    th=-1;      % prevent any stopping if l has no fractional part
end
gg=zeros(l+1,1);
im=repmat(1:k,1,nb); im=im(:);

% index arrays for replication

wk=ones(k,1);
wp=ones(1,p);
% wn=ones(1,n);
%
% % Main calculation loop
%
% We have the following relationships to [1] where i and k index
% the data values and cluster centres respectively:
%
%   This program     [Zhang2000]                            Equation  
%
%     d(i,:)            x_i                                 input data
%     x(k,:)            m_k                                 cluster centres
%     py(k,i)           (d_ik)^2
%     dm(i)'            d_i,min^2
%     pr(k,i)           (d_i,min/d_ik)^2
%     pe(k,i)           (d_i,min/d_ik)^p                    (7.6) 
%     qik(k,i)          q_ik                                (7.2)
%     qk(k)             q_k                                 (7.3)
%     qik(k,i)./qk(k)   p_ik                                (7.4)
%     se(i)'            d_i,min^p * sumk(d_ik^-p)
%     xf(i)'            d_i,min^-2 / sumk(d_ik^-p)
%     xg(i)'            d_i,min^-(p+2) / sumk(d_ik^-p)^2


ss=sd+1;        % one extra loop at the start
g=0;                % dummy initial value of g
xn=zeros(n,1);
for j=1:l

    g1=g;                           % save old performance
    x1=x;                           % save old centres
    % first do partial chunk

    jx=n-(nl-1)*nb;
    ii=1:jx;
    kx=repmat(ii,k,1);
    km=repmat(1:k,1,jx);
    py=reshape(sum((d(kx(:),:)-x(km(:),:)).^2,2),k,jx);
    [dm,xn(ii)]=min(py,[],1);                 % min value in each column gives nearest centre
    dmk=dm(wk,:);                   % expand into a matrix
    dq=py>dmk;                      % update only these values
    pr=ones(k,jx);                   % leaving others at 1
    pr(dq)=dmk(dq)./py(dq);            % ratio of min(py)./py
    pe=pr.^eh;
    se=sum(pe,1);
    xf=dm.^(eh-1)./se;
    g=xf*dm.';                     % performance criterion (divided by k)
    xg=xf./se;
    qik=xg(wk,:).*pe.*pr;           % qik(k,i) is equal to q_ik in [Zhang2000]
    qk=sum(qik,2);
    xs=qik*d(ii,:);
    ix=jx+1;
    for il=2:nl
        jx=jx+nb;        % increment upper limit
        ii=ix:jx;
        kx=ii(wk,:);
        py=reshape(sum((d(kx(:),:)-x(im,:)).^2,2),k,nb);
        [dm,xn(ii)]=min(py,[],1);                 % min value in each column gives nearest centre
        dmk=dm(wk,:);                   % expand into a matrix
        dq=py>dmk;                      % update only these values
        pr=ones(k,nb);                   % leaving others at 1
        pr(dq)=dmk(dq)./py(dq);            % ratio of min(py)./py
        pe=pr.^eh;
        se=sum(pe,1);
        xf=dm.^(eh-1)./se;
        g=g+xf*dm.';                     % performance criterion (divided by k)
        xg=xf./se;
        qik=xg(wk,:).*pe.*pr;           % qik(k,i) is equal to q_ik in [Zhang2000]
        qk=qk+sum(qik,2);
        xs=xs+qik*d(ii,:);
        ix=jx+1;
    end
    gg(j)=g;
    x=xs./qk(:,wp);
    if g1-g<=th*g1
        ss=ss-1;
        if ~ss break; end  %  stop if improvement < threshold for sd consecutive iterations
    else
        ss=sd;
    end
end
gg=gg(1:j)*k/n;                       % scale and trim the performance criterion vector
g=g(end);
% gg' % *** DEBUIG ***
if nargout>1
    x=x1;                               % go back to the previous x values if G and/or XN value is output
end

