function [mz,vz,wz]=gaussmixd(y,m,v,w,a,b,f,g)
%GAUSSMIXD marginal and conditional Gaussian mixture densities
%
% Usage: (1) [mz,vz,wz]=gaussmixd(y,m,v,w); % If {m,v,w} has dimension p, then y specifies
%                                           % elements 1:q and {mz,vz,wz} is a GMM for elements q+1:p
%        (2) [mz,vz,wz]=gaussmixd(y,m,v,w,[],[1 3 4]);	% y specifies elements 1,3,4 and {mz,vz,wz}
%                                                       % is a GMM for the remaining elements
%
% Inputs: Input mixture: k mixtures, p dimensions
%         Output mixture: k mixtures, r dimensions
%         Conditioning data: n data values, q dimensions
%
%   Y(n,q) = conditioning input data: x*a'+ b'= y
%   M(k,p) = mixture means for x(p)
%   V(k,p) or V(p,p,k) variances (diagonal or full)
%   W(k,1) = mixture weights
%   A(q,p) = conditioning transformation: y=x*A'+ B' (where y and x are row vectors).
%   B(q,1)   If A is omitted or null, y=x*I(B,:)' where I is the identity matrix.
%            If B is also omitted or null, y=x*I(1:q,:)'.
%   F(r,p) = output transformation z = x*F'+G'  (where z and x are row vectors).
%   G(r,1)   If F is omitted or null, z = x*f' where I is the identity matrix.
%            If G is also omitted or null, z=x*I(q+1:p,:)' or, if A was also null,
%            the complement of y.
%
% Outputs (if 1 or 2 outputs specified):
%
%   MZ(n,r) = Global mean of z for each y
%   VZ(r,r,n) Global full covariances of z (now dependent on y)
%
% Outputs (if 3 outputs specified):
%
%   MZ(k,r,n) = mixture means of z for each y
%   VZ(k,r) or VZ(r,r,k) variances of z (diagonal or full); surprisingly it is independent of y
%   WZ(n,k)
%
% The output mixture covariances will be diagonal if either r=1 or else the following three
% conditions are all true:
%  (a) the input mixture covariances are diagonal and
%  (b) matrix A has at most one non-zero element in any row or column and
%  (c) matrix F has at most one non-zero element in any column
%
% Several of the output variables can be squeezed if r=1 but this is not done automatically.
%
% This routine can be used for inference: if you train a GMM on p-dimensional data
% then, if y(n,q) contains observations of the first q components, then z=gaussmixd(y,m,v,w)
% will return the estimated values of the remaining p-q components.
%
% See also: gaussmix, gaussmixg, gaussmixp, randvec

%      Copyright (C) Mike Brookes 2000-2012
%      Version: $Id: gaussmixd.m 7784 2016-04-15 11:09:50Z dmb $
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

[n,q]=size(y);
[k,p]=size(m);
% We set fv=1 if input V contains full covariance matrices with dimensions V(r,r,k).
% This is true either if V has >2 dimensions or V=V(r,r) and r>k
fv=ndims(v)>2 || size(v,1)>k; % full covariance matrix is supplied
anull=nargin<5 || isempty(a); % no A matrix specified
if anull
    a=eye(p);
    if nargin<6 || isempty(b)
        b=1:q;
    end
    a=a(b,:);
    b0=b;  % save row selection
    b=zeros(q,1);
elseif nargin<6 || isempty(b)
    b=zeros(q,1);
end
if nargin<7 || isempty(f)
    f=eye(p);
    if nargin<8 || isempty(g)
        if anull
            f(b0,:)=[];
        else
            f=f(q+1:end,:);
        end
    else
        f=f(g,:);  % G selects the output variables
    end
    r=size(f,1);
    g=zeros(r,1);
elseif nargin<8 || isempty(g)
    r=size(f,1);
    g=zeros(r,1);
else
    r=size(f,1);
end

yb=y-repmat(b',n,1);                 % remove the b term
[lp,wz]=gaussmixp(yb,m,v,w,a);     % find mixture weights
mz=zeros(n,r,k);
ma=a~=0;                            % check for sparse a and f matrices (common case)
mf=f~=0;
% We set dvo=1 if the output mixture covariances are structurally
% diagonal. This is the case if either:
% (1) r=1 since in this case they are scalar values, or else
% (2) the following three conditions are all true:
%     (a) the input mixture covariances are diagonal and
%     (b) matrix A has at most one non-zero element in any row or column and
%     (c) matrix F has at most one non-zero element in any column
dvo=r==1 || (~fv && all(sum(ma,1)<=1) && all(sum(ma,2)<=1) && all(sum(mf,1)<=1));
if dvo    % structurally diagonal output covariances
    vz=zeros(k,r); % diagonal output variances (one row per mixture) independent of y
    for i=1:k                               % loop for each mixture
        if fv
            vi=v(:,:,i);
        else
            vi=diag(v(i,:));                % convert to full covariance matrices
        end
        hi=vi*a'/(a*vi*a');                 % regression coefficient matrix (p#q]
        vz(i,:)=diag(f*(vi-hi*a*vi)*f')';   % variance of z (independent of y)
        mi=m(i,:);                          % input mean for mixure i
        m0=(mi-mi*a'*hi')*f'+g';            % y-independent part of mean
        mz(:,:,i)=(repmat(m0,n,1)+yb*hi'*f'); % mean for each y
    end
else
    vz=zeros(r,r,k);
    for i=1:k                               % loop for each mixture
        if fv
            vi=v(:,:,i);
        else
            vi=diag(v(i,:));                % convert to full covariance matrices
        end
        hi=vi*a'/(a*vi*a');                 % regression coefficient matrix (p#q]
        vz(:,:,i)=f*(vi-hi*a*vi)*f';        % variance of z (independent of y)
        mi=m(i,:);                          % input mean for mixure i
        m0=(mi-mi*a'*hi')*f'+g';            % y-independent part of mean
        mz(:,:,i)=(repmat(m0,n,1)+yb*hi'*f'); % mean for each y
    end
end
if nargout<3
    mt=reshape(sum(reshape(mz,n*r,k).*repmat(wz,r,1),2),n,r);       % global mean
    if nargout>1        % need to calculate global variance as well
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % We calculate the global covariance by adding up the weighted mixture
        % covariances corrected for the fact that the mixture mean does not equal
        % the global mean.
        % To save calculations, we calculate only the lower triangle of the
        % symmetric covariance matrix and then expand it at the end into a full matrix.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rl=r*(r+1)/2;   % number of elements in lower triangular covariance matrix
        lix=1:r^2;
        cix=repmat(1:r,r,1);
        rix=cix';
        lix(cix>rix)=[];                                        % index of lower triangular elements
        rlix=rix(lix);
        clix=cix(lix);
        lixi=zeros(r,r);
        lixi(lix)=1:rl;
        lixi=lixi';
        lixi(lix)=1:rl;                                        % reverse index to build full matrices

        vt=zeros(n,rl);                % reserve space for lower triangular output covariances
        for i=1:k
            if dvo
                vi=diag(vz(i,:));
            else
                vi=vz(:,:,i);
            end
            mzt=mz(:,:,i)-mt;
            vt=vt+repmat(wz(:,i),1,r^2).*(repmat(vi(lix),n,1)+mzt(:,rlix).*mzt(:,clix));
        end
        vz=permute(reshape(vt(:,lixi),[n,r,r]),[2 3 1]);
    end
    mz=mt;
else
    mz=permute(mz,3:-1:1);
end
