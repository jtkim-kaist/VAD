function [u,v,p,r] = maxgauss(m,c,d)
%MAXGAUSS determine gaussian approximation to max of a gaussian vector [p,u,v,r]=(m,c,d)
%
% Inputs:
%       m(N,1) is the mean vector of length N
%       c(N,N) is Covariance matrix (or c(N,1) is vector of covariances)
%       d(N,K) is Covariance w.r.t some other variables
%
% Outputs:
%       u is mean(max(x)) where x is the random variable of length N
%       v is var(max(x))
%       p(N,1) is prob of each element being the max
%       r(1,K) is covariance between max(x) and the variables corresponding to the columns of d
%
% To find the min instead of max, just negate m and u.

% The algorithm combines the elements of m in pairs in sequence as suggested
% in [1]. Errors are introduced because max(x(i),x(j)) is wrongly assumed to be
% gaussian. To minimize the errors, we use a greedy algorithm (approximately as
% in [2]) in which the chosen pair has the greatest imbalance in selection probability.
%
%
% Refs: [1] C. E. Clark. The greatest of a finite set of random variables.
%           Operations Research, 9(2):145–162, March 1961.
%       [2] D. Sinha, H. Zhou, and N. V. Shenoy.
%           Advances in Computation of the Maximum of a Set of Gaussian Random Variables.
%           IEEE Trans Computer-Aided Design of Integrated Circuits and Systems,
%           26(8):1522–1533, Aug. 2007. doi: 10.1109/TCAD.2007.893544.

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: maxgauss.m 713 2011-10-16 14:45:43Z dmb $
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

nrh=-sqrt(0.5);
kpd=sqrt(0.5/pi);
m=m(:);
nm=length(m);
if nargin<2
    c=eye(nm);
elseif numel(c)==nm
    c=diag(c);
end
p=eye(nm);
ix=find(m==-Inf);       % remove negative infinities immediately
if ~isempty(ix)
    m(ix)=[];
    c(ix,:)=[];
    c(:,ix)=[];
    p(:,ix)=[];
    nm=length(m);
end
ix=0:nm^2-1;                        % index for nm x nm matrix elements (from 0)
ix = ix(ix<(nm+1)*floor(ix/nm));	% only keep the strict upper triangle
mij=zeros(2,1);                     % store for means
while nm>1

    % calculate scaled differences in means

    gm=m(:,ones(1,nm));
    gm = gm-gm.';               % find the difference between all pairs of means
    gm=gm(ix+1);                % only keep the strict upper triangular elements
    cd=diag(c);
    gv=cd(:,ones(1,nm))-c;
    gv=gv+gv.';
    gv=gv(ix+1);                % These are the corresponding variance sums
    jx=find(gv<=0);
    if ~isempty(jx)             % special case: two variables differ by a constant
        jx=jx(1);               % take first pair for which this is true
        j=floor(ix(jx)/nm);
        i=ix(jx)-nm*j+1;
        j=j+1;
        dm=gm(jx);
        if dm>0                 % if x(i)>x(j) then abolish j
            m(j)=[];
            c(j,:)=[];
            c(:,j)=[];
            p(:,j)=[];
        else                    % if x(i)<=x(j) then abolish i
            m(i)=[];
            c(i,:)=[];
            c(:,i)=[];
            p(:,i)=[];
        end
    else
        % select the pair of variables with the the highest ratio of
        % squared difference in means to sum of variances and combine into a single variable
        [gg,jx]=max(gm.^2./gv);         % jx indicates which pair to combine
        j=floor(ix(jx)/nm);             % convert jx into (i,j) pair
        i=ix(jx)-nm*j+1;
        j=j+1;                          % combine variables i and j
        dm=gm(jx);                      % mean of x(i)-x(j)
        ds=sqrt(gv(jx));                % std dev of x(i)-x(j)
        dms=dm/ds;
        q = 0.5 * erfc(nrh*dms);        % q =normcdf(dm,0,ds) = Phi(dms) = prob x(i) > x(j)
        f=kpd*exp(-0.5*dms^2);          % f=phi(dms)=ds*normpdf(dm,0,ds)
        mij(1)=m(i);
        mij(2)=m(j);
        u=dm*q+mij(2)+ds*f;                                         % mean of max{x(i),x(j)}
        v=(mij(1)+mij(2)-u)*u +cd(i)*q+cd(j)*(1-q)-mij(1)*mij(2);   % variance of max{x(i),x(j)}

        % replace x(i) with max{x(i),x(j)}

        m(i)=u;
        c(i,:)=q*c(i,:)+(1-q)*c(j,:);
        c(:,i)=c(i,:).';
        c(i,i)=v;
        p(:,i)=q*p(:,i)+(1-q)*p(:,j);
        
        % now abolish x(j)
        
        m(j)=[];
        c(j,:)=[];
        c(:,j)=[];
        p(:,j)=[];
    end

    ix=ix(1:(nm-1)*(nm-2)/2);
    ix=ix-floor(ix/nm);
    nm=nm-1;
end         % main while loop
u=m(1);
v=c(1);
p=p/sum(p);     % force sum=1
if nargin>2
    r=p.'*d;
end