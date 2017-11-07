function [mm,mc]=gaussmixm(m,v,w,z)
% GAUSSMIXM estimate mean and variance of the magnitude of a GMM
%
%  Inputs:  M(K,P)   is the mean vectors (one row per mixture)
%           V(K,P)   are diagonal covariances (one row per mixture) [ones(K,P)]
%        or V(P,P,K) are full covariance matrices (one per mixture)
%           W(K)     are the mixture weights [ones(K,1)/K]
%           Z(T,P)   each row is an origin to measure magnitude from [zeros(1,P)]
%
% Outputs:  MM(T,1)  mean of sqrt((x-z(t))'*(x-z(t))) where x is the GMM random variate
%           MC(T,T)  covariance matrix of sqrt((x-z(t))'*(x-z(t)))
%
% We approximate the magnitude of each mixture as a Nakagami-m distribution and
% match the moments of its squared variate. We approximate the normalized 
% correlation matrix of |x| (i.e. with unit diagonal) to be the same as that
% of |x|^2 which we can calculate exactly.
% Answers are exact for P=1 or when all M=0. Accuracy improves otherwise for
% large P or M.

%      Copyright (C) Mike Brookes 2015
%      Version: $Id: gaussmixm.m 5957 2015-03-26 16:59:23Z dmb $
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
[k,p]=size(m);                          % k = # mixtures, p = vector dimension
if nargin<4 || ~numel(z)
    z=zeros(1,p);
end
if nargin<3 || ~numel(w)
    w=ones(k,1);
end% default to uniform weights
if nargin<2 || ~numel(v)
    v=ones(k,p);                    % default to unit variance
end
[t,p1]=size(z);
if p~=p1 || (p>2 && t>1 && nargout>2)
    error('unable to calculate a covariance matrix');
end
w=w(:)/sum(w);                          % make the weights sum to 1 and force column vector
if p==1                                 % treat 1D case specially since we have an exact formula
    s=sqrt(v(:));                       % normalize mixture means by the std dev
    mt=m(:,ones(1,t))-z(:,ones(1,k))';  % shifted mean for each mixture x origin
    mts=mt./s(:,ones(1,t));             % shifted mean for each mixture x origin normalized to unit SD
    ncdf=normcdf(-mts);
    npdf=normpdf(-mts);
    mm=(mts.*(1-2*ncdf)+2*npdf)'*(s.*w); % exact mean of |X|
    mc=diag((mts.^2+1)'*(v(:).*w)); % diagonal variance elements in case t==1
    if nargout>1 && t>1                 % we need to calculate a covariance matrix
        for it=1:t
            for jt=1:it-1
                mc(it,jt)=w.'*((v(:)+mt(:,it).*mt(:,jt)).*(1-2*abs(ncdf(:,it)-ncdf(:,jt)))+ ...
                    2*s.*sign(mt(:,jt)-mt(:,it)).*(npdf(:,it).*mt(:,jt)-npdf(:,jt).*mt(:,it)));
                mc(jt,it)=mc(it,jt);
            end
        end   
    end
    mc=mc-mm*mm';                       % convert to variance
else                                        % p>1 case
    fullv=ndims(v)>2 || size(v,1)>k;      	% full covariance matrix is supplied
    if fullv                                % full covariance matrix is supplied
        ms=repmat(sum(m.^2+v(repmat((1:p+1:p^2)',1,k)+repmat(0:p^2:(k-1)*p^2,p,1))',2),1,t)-2*m*z'+repmat(sum(z.^2,2)',k,1);  % mean of |x|^2 for each mixture x origin
        vsf=zeros(t,t,k);
        vs=zeros(k,t);
        for i=1:k
            zmi=repmat(m(i,:),t,1)-z;
            si=v(:,:,i);
            vsi=2*trace(si^2)+4*zmi*si*zmi';
            vsf(:,:,i)=vsi;
            vs(i,:)=diag(vsi);
        end
    else                                        % else diagonal covariance matrix supplied
        ms=repmat(sum(m.^2,2)+sum(v,2),1,t)-2*m*z'+repmat(sum(z.^2,2)',k,1);  % mean of |x|^2 for each mixture x origin
        vsc=sum(v.*(2*v+4*m.^2),2);             % origin-independent part of  |x|^2 variance for each mixture
        vmz=4*(v.*m)*z';                        % origin-linear part of  |x|^2 variannce for each mixture x origin
        vs=repmat(vsc,1,t)-2*vmz+4*v*z.^2';     % variance of |x|^2 for each mixture x origin
    end
    nm=ms.^2./vs;                               % Nakagami-m parameter per mixture x origin
    mmk=(exp(gammaln(nm+0.5)-gammaln(nm)).*sqrt(ms./nm)); % mean of Nakagami-m distrbution per mixture x origin
    mm=mmk'*w;                                  % mean per origin
    mc = ms'*w-mm.^2 ;                          % variance in case t==1
    if nargout>1 && t>1                         % we need to calculate the t x t covariance matrix w.r.t. origins
        mvksd=sqrt(ms-mmk.^2);                  % std deviation of |x| per mixture x origin
        mc=zeros(t,t);
        if fullv                                % full covariance matrix is supplied
            for i=1:k
                vsisd=sqrt(vs(i,:));
                mcm=(mvksd(i,:)'*mvksd(i,:)).*vsf(:,:,i)./(vsisd'*vsisd)+ mmk(i,:)'*mmk(i,:); % estimated |x| correlation matrix for mixture i
                mc=mc+w(i)*mcm;             % add onto total correlation matrix
            end
        else                                % else diagonal covariance matrix supplied
            for i=1:k                       % loop for each mixture component
                iit=repmat(i,t,1);
                vsi=vsc(i)-vmz(iit,:)-vmz(iit,:)'+4*z.*v(iit,:)*z'; % covariance matriz for |x|^2 for mixture i
                vsisd=sqrt(vs(i,:));
                mcm=(mvksd(i,:)'*mvksd(i,:)).*vsi./(vsisd'*vsisd)+ mmk(i,:)'*mmk(i,:); % estimated |x| correlation matrix for mixture i
                mc=mc+w(i)*mcm;   % add onto total correlation matrix
            end
        end
        mc=mc-mm*mm'; % convert to covariance matrix
    end
end