function x=randvec(n,m,c,w,mode)
%RANDVEC  Generate real or complex GMM/lognormal random vectors X=(N,M,C,W,MODE)
% generates a random matrix of size (|n|,p) where p is the maximum
% dimension of M or C
%  Inputs:  N        is the number of points to generate
%           M(K,P)   is the mean vectors (one row per mixture)
%           C(K,P)   are diagonal covariances (one row per mixture)
%        or C(P,P,K) are full covariance matrices (one per mixture)
%           W(K)     are the mixture weights (or omit if all mixtures have equal weight)
%           MODE     character string specifying options:
%                       g = real gaussian [default]
%                       c = complex gaussian
%                       l = lognormal
%
% Outputs:  X(N,P) is the output data
%
% Note that cov(x) = E(x'*x) - m'*m where x and m are row vectors (for complex data this
% may be the conjugate of what you expect but agrees with MATLAB's COV() function.

% Bugs/suggestions
%    (1)  New mode 'x' to approximate chi squared

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: randvec.m 5453 2014-11-19 13:10:51Z dmb $
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

% first sort out the input arguments

sm=size(m);
if nargin<3
    c=ones(sm);         % default to unit variance
end
sc=size(c);
p=max(sm(2),sc(2));     % data dimension
k=sm(1);                % number of mixtures
fullc=(length(sc)>2) || (sc(1)>k);
if nargin<4
    mode='g';   % default to gaussian
    w=ones(k,1);
else
    if ischar(w)
        mode = w;       % w argument has been omitted
        w=ones(k,1);
    elseif nargin<5
        mode='g';
    end
end
ty=mode(1);   % ignore all but first character for type
x=zeros(n,p);   % initialize output array
if sm(2)~=p
    m=repmat(m,1,p);    % if m is given as a scalar
end
if sc(2) ~=p
    c=repmat(c,1,p);    % if c is given as a scalar
end
q=sqrt(0.5);
if k>1
    kx=randiscr(w,n);
else
    kx=ones(n,1);
end
for kk=1:k
    nx=find(kx==kk);
    nn=length(nx);
    if nn       % check if we need to generate any from mixture kk

        % extract the mean and cov for this mixture

        mm=m(kk,:);     % mean vector
        if fullc        % full covariance matrix
            cc=c(:,:,kk);
            if ty=='l'      % lognormal distribution - convert mean and covariance
                cc=log(1+cc./(mm.'*mm));
                mm=log(mm)-0.5*diag(cc).';
            end
        else
            cc=c(kk,:);
            if ty=='l'      % lognormal distribution - convert mean and covariance
                cc=log(1+cc(:).'./mm.^2);
                mm=log(mm)-0.5*cc;
            end
        end

        % now generate nn complex or real values

        if ty=='c'  % generate complex values
            xx=q*randn(nn,p)+1i*q*randn(nn,p); % complex-valued unit variance values
        else
            xx=randn(nn,p);   % real-valued unit variance values
        end;

        % scale by the square root of the covariance matrix

        if fullc   % full covariance covariance
            [v,d]=eig((cc+cc')/2);   % force covariance matrix to be hermitian
            xx=(xx.*repmat(sqrt(abs(diag(d))).',nn,1))*v'+repmat(mm,nn,1); % and also positive definite
        else
            xx=xx.*repmat(sqrt(abs(cc)),nn,1)+repmat(mm,nn,1); % different mean/cov for each column
        end
        x(nx,:)=xx;
    end
end
if ty=='l'  % lognormal distribution
    x=exp(x);
end
if ~nargout
    if p==1
        if ty=='c'
            plot(real(x), imag(x),'+');
            xlabel('Real');
            ylabel('Imag');
        else
            nbin=max(min(floor(sqrt(n)),50),5);
            hist(x,nbin);
            xlabel('X');
            ylabel('Frequency');
        end
    else
        [vv,iv]=sort(var(x,0,1));
        iv=sort(iv(end-1:end));
        plot(real(x(:,iv(1))), real(x(:,iv(2))),'+');
        if ty=='c'
            xylab='Real[ x(%d) ]';
        else
            xylab='x(%d)';
        end
        xlabel(sprintf(xylab,iv(1)));
        ylabel(sprintf(xylab,iv(2)));
    end
end

