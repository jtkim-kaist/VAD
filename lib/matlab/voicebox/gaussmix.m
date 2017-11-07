function [m,v,w,g,f,pp,gg]=gaussmix(x,c,l,m0,v0,w0,wx)
%GAUSSMIX fits a gaussian mixture pdf to a set of data observations [m,v,w,g,f]=(x,c,l,m0,v0,w0,wx)
%
% Usage:
%    (1) [m,v,w]=gaussmix(x,[],[],k);    % create GMM with k mixtures and diagonal covariances
%    (2) [m,v,w]=gaussmix(x,[],[],k,'v);    % create GMM with k mixtures and full covariances
%
% Inputs: n data values, k mixtures, p parameters, l loops
%
%     X(n,p)   Input data vectors, one per row.
%     c(1)     Minimum variance of normalized data (Use [] to take default value of 1/n^2)
%     L        The integer portion of l gives a maximum loop count. The fractional portion gives
%              an optional stopping threshold. Iteration will cease if the increase in
%              log likelihood density per data point is less than this value. Thus l=10.001 will
%              stop after 10 iterations or when the increase in log likelihood falls below
%              0.001.
%              As a special case, if L=0, then the first three outputs are omitted.
%              Use [] to take default value of 100.0001
%     M0       Number of mixtures required (or initial mixture means - see below)
%     V0       Initialization mode:
%                'f'    Initialize with K randomly selected data points [default]
%                'p'    Initialize with centroids and variances of random partitions
%                'k'    k-means algorithm ('kf' and 'kp' determine initialization)
%                'h'    k-harmonic means algorithm ('hf' and 'hp' determine initialization) [default]
%                's'    do not scale data during initialization to have equal variances
%                'm'    M0 contains the initial centres
%                'v'    full covariance matrices
%              Mode 'hf' [the default] generally gives the best results but 'f' is faster and often OK
%     W0(n,1)  Data point weights
%
%     Alternatively, initial values for M0, V0 and W0 can be given  explicitly:
%
%     M0(k,p)  Initial mixture means, one row per mixture.
%     V0(k,p)  Initial mixture variances, one row per mixture.
%      or V0(p,p,k)  one full-covariance matrix per mixture
%     W0(k,1)  Initial mixture weights, one per mixture. The weights should sum to unity.
%     WX(n,1)  Data point weights
%
% Outputs: (Note that M, V and W are omitted if L==0)
%
%     M(k,p)   Mixture means, one row per mixture. (omitted if L==0)
%     V(k,p)   Mixture variances, one row per mixture. (omitted if L==0)
%       or V(p,p,k) if full covariance matrices in use (i.e. either 'v' option or V0(p,p,k) specified)
%     W(k,1)   Mixture weights, one per mixture. The weights will sum to unity. (omitted if L==0)
%     G       Average log probability of the input data points.
%     F        Fisher's Discriminant measures how well the data divides into classes.
%              It is the ratio of the between-mixture variance to the average mixture variance: a
%              high value means the classes (mixtures) are well separated.
%     PP(n,1)  Log probability of each data point
%     GG(l+1,1) Average log probabilities at the beginning of each iteration and at the end
%
% The fitting procedure uses one of several initialization methods to create an initial guess
% for the mixture centres and then uses the EM (expectation-maximization) algorithm to refine
% the guess. Although the EM algorithm is deterministic, the initialization procedures use 
% random numbers and so the routine will not give identical answers if you call it multiple
% times with the same input data.

%  Bugs/Suggestions
%     (1) Allow processing in chunks by outputting/reinputting an array of sufficient statistics
%     (2) Other initialization options:
%              'l'    LBG algorithm
%              'm'    Move-means (dog-rabbit) algorithm
%     (3) Allow updating of weights-only, not means/variances

%      Copyright (C) Mike Brookes 2000-2009
%      Version: $Id: gaussmix.m 7784 2016-04-15 11:09:50Z dmb $
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

[n,p]=size(x);
wn=ones(n,1);
mx0=sum(x,1)/n;         % calculate mean and variance of input data in each dimension
vx0=sum(x.^2,1)/n-mx0.^2;
sx0=sqrt(vx0);
sx0(sx0==0)=1;      % do not divide by zero when scaling
scaled=0;           % data is not yet scaled
memsize=voicebox('memsize');    % set memory size to use
if isempty(c)
    c=1/n^2;
else
    c=c(1);         % just to prevent legacy code failing
end
fulliv=0;           % initial variance is not full
if isempty(l)
    l=100+1e-4;         % max loop count + stopping threshold
end
if nargin<5 || isempty(v0) || ischar(v0)             % no initial values specified for m0, v0, w0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  No initialvalues given, so we must use k-means or equivalent
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<6
        if nargin<5 || isempty(v0)
            v0='hf';                 % default initialization mode: hf
        end
        wx=wn;                      % no data point weights
    else
        wx=w0(:);                   % data point weights
    end
    if any(v0=='m')
        k=size(m0,1);
    else
        k=m0;
    end
    fv=any(v0=='v');                % full covariance matrices requested
    if n<=k                         % each data point can have its own mixture
        xs=(x-mx0(wn,:))./sx0(wn,:);          % scale the data
        m=xs(mod((1:k)-1,n)+1,:);   % just include all points several times
        v=zeros(k,p);               % will be set to floor later
        w=zeros(k,1);
        w(1:n)=1/n;
        if l>0
            l=0.1;                  % no point in iterating
        end
    else                            % more points than mixtures
        if any(v0=='s')
            xs=x;                   % do not scale data during initialization
        else
            xs=(x-mx0(wn,:))./sx0(wn,:);  % else scale now
            if any(v0=='m')
                m=(m0-mx0(ones(k,1),:))./sx0(ones(k,1),:);  % scale specified means as well
            end
        end
        w=repmat(1/k,k,1);                      % all mixtures equally likely
        if any(v0=='k')                         % k-means initialization
            if any(v0=='m')
                [m,e,j]=v_kmeans(xs,k,m);
            elseif any(v0=='p')
                [m,e,j]=v_kmeans(xs,k,'p');
            else
                [m,e,j]=v_kmeans(xs,k,'f');
            end
        elseif any(v0=='h')                     % k-harmonic means initialization
            if any(v0=='m')
                [m,e,j]=kmeanhar(xs,k,[],4,m);
            else
                if any(v0=='p')
                    [m,e,j]=kmeanhar(xs,k,[],4,'p');
                else
                    [m,e,j]=kmeanhar(xs,k,[],4,'f');
                end
            end
        elseif any(v0=='p')                     % Initialize using a random partition
            j=ceil(rand(n,1)*k);                % allocate to random clusters
            j(rnsubset(k,n))=1:k;               % but force at least one point per cluster
            for i=1:k
                m(i,:)=mean(xs(j==i,:),1);
            end
        else
            if any(v0=='m')
                m=m0;                           % use specified centres
            else
                m=xs(rnsubset(k,n),:);          % Forgy initialization: sample k centres without replacement [default]
            end
            [e,j]=v_kmeans(xs,k,m,0);             % find out the cluster allocation
        end
        if any(v0=='s')
            xs=(x-mx0(wn,:))./sx0(wn,:);      % scale data now if not done previously
        end
        v=zeros(k,p);                   % diagonal covariances
        w=zeros(k,1);
        for i=1:k
            ni=sum(j==i);               % number assigned to this centre
            w(i)=(ni+1)/(n+k);          % weight of this mixture
            if ni
                v(i,:)=sum((xs(j==i,:)-repmat(m(i,:),ni,1)).^2,1)/ni;
            else
                v(i,:)=zeros(1,p);
            end
        end
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%
    % use initial values given as input parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [k,p]=size(m0);
    xs=(x-mx0(wn,:))./sx0(wn,:);          % scale the data
    m=(m0-mx0(ones(k,1),:))./sx0(ones(k,1),:);          % and the means
    v=v0;
    w=w0;
    fv=ndims(v)>2 || size(v,1)>k;                       % full covariance matrix is supplied
    if fv
        mk=eye(p)==0;                                    % off-diagonal elements
        fulliv=any(v(repmat(mk,[1 1 k]))~=0);            % check if any are non-zero
        if ~fulliv
            v=reshape(v(repmat(~mk,[1 1 k])),p,k)'./repmat(sx0.^2,k,1);   % just pick out and scale the diagonal elements for now
        else
            v=v./repmat(sx0'*sx0,[1 1 k]);              % scale the full covariance matrix
        end
    end
    if nargin<7
        wx=wn;              % no data point weights
    end
end
if length(wx)~=n
    error('%d datapoints but %d weights',n,length(wx));
end
lsx=sum(log(sx0));
xsw=xs.*repmat(wx,1,p); % weighted data points
nwt=sum(wx);        % number of data points counting duplicates
if ~fulliv          % initializing with diagonal covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagonal Covariance matrices  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v=max(v,c);         % apply the lower bound
    xs2=xs.^2.*repmat(wx,1,p);          % square and weight the data for variance calculations

    % If data size is large then do calculations in chunks

    nb=min(n,max(1,floor(memsize/(8*p*k))));    % chunk size for testing data points
    nl=ceil(n/nb);                  % number of chunks
    jx0=n-(nl-1)*nb;                % size of first chunk

    im=repmat(1:k,1,nb); im=im(:);
    th=(l-floor(l))*n;
    sd=(nargout > 3*(l~=0)); % = 1 if we are outputting log likelihood values
    lp=floor(l)+sd;   % extra loop needed to calculate final G value

    lpx=zeros(1,n);             % log probability of each data point
    wk=ones(k,1);
    wp=ones(1,p);
    wnb=ones(1,nb);
    wnj=ones(1,jx0);

    % EM loop

    g=0;                           % dummy initial value for comparison
    gg=zeros(lp+1,1);
    ss=sd;                       % initialize stopping count (0 or 1)
    for j=1:lp
        g1=g;                    % save previous log likelihood (2*pi factor omitted)
        m1=m;                       % save previous means, variances and weights
        v1=v;
        w1=w;
        vi=-0.5*v.^(-1);                % data-independent scale factor in exponent
        lvm=log(w)-0.5*sum(log(v),2);   % log of external scale factor (excluding -0.5*p*log(2pi) term)

        % first do partial chunk (of length jx0)

        jx=jx0;
        ii=1:jx;                        % indices of data points in this chunk
        kk=repmat(ii,k,1);              % kk(jx,k): one row per data point, one column per mixture
        km=repmat(1:k,1,jx);            % km(jx,k): one row per data point, one column per mixture
        py=reshape(sum((xs(kk(:),:)-m(km(:),:)).^2.*vi(km(:),:),2),k,jx)+lvm(:,wnj); % py(k,jx) pdf of each point with each mixture
        mx=max(py,[],1);                % mx(1,jx) find normalizing factor for each data point to prevent underflow when using exp()
        px=exp(py-mx(wk,:));            % find normalized probability of each mixture for each datapoint
        ps=sum(px,1);                   % total normalized likelihood of each data point
        px=px./ps(wk,:);                % relative mixture probabilities for each data point (columns sum to 1)
        lpx(ii)=log(ps)+mx;
        pk=px*wx(ii);                       % pk(k,1) effective number of data points for each mixture (could be zero due to underflow)
        sx=px*xsw(ii,:);
        sx2=px*xs2(ii,:);
        for il=2:nl                     % process the data points in chunks
            ix=jx+1;
            jx=jx+nb;                   % increment upper limit
            ii=ix:jx;                   % indices of data points in this chunk
            kk=repmat(ii,k,1);
            py=reshape(sum((xs(kk(:),:)-m(im,:)).^2.*vi(im,:),2),k,nb)+lvm(:,wnb);
            mx=max(py,[],1);            % find normalizing factor for each data point to prevent underflow when using exp()
            px=exp(py-mx(wk,:));        % find normalized probability of each mixture for each datapoint
            ps=sum(px,1);               % total normalized likelihood of each data point
            px=px./ps(wk,:);            % relative mixture probabilities for each data point (columns sum to 1)
            lpx(ii)=log(ps)+mx;
            pk=pk+px*wx(ii);                % pk(k,1) effective number of data points for each mixture (could be zero due to underflow)
            sx=sx+px*xsw(ii,:);
            sx2=sx2+px*xs2(ii,:);
        end
        g=lpx*wx;                       % total log probability summed over all data points
        gg(j)=g;                        % save log prob at each iteration
        w=pk/nwt;                       % normalize to get the weights
        if pk                           % if all elements of pk are non-zero
            m=sx./pk(:,wp);             % calculate mixture means
            v=sx2./pk(:,wp);            % and variances
        else
            wm=pk==0;                   % mask indicating mixtures with zero weights
            nz=sum(wm);              	% number of zero-weight mixtures
            [vv,mk]=sort(lpx);          % find the lowest probability data points
            m=zeros(k,p);               % initialize means and variances to zero (variances are floored later)
            v=m;
            m(wm,:)=xs(mk(1:nz),:); 	% set zero-weight mixture means to worst-fitted data points
            w(wm)=1/n;               	% set these weights non-zero
            w=w*n/(n+nz);            	% normalize so the weights sum to unity
            wm=~wm;                 	% mask for non-zero weights
            m(wm,:)=sx(wm,:)./pk(wm,wp);  % recalculate means and variances for mixtures with a non-zero weight
            v(wm,:)=sx2(wm,:)./pk(wm,wp);
        end
        v=max(v-m.^2,c);                % apply floor to variances
        if g-g1<=th && j>1
            if ~ss, break; end  %  stop
            ss=ss-1;       % stop next time
        end

    end
    if sd && ~fv  % we need to calculate the final probabilities
        pp=lpx'-0.5*p*log(2*pi)-lsx;   % log of total probability of each data point
        gg=gg(1:j)/n-0.5*p*log(2*pi)-lsx;    % average log prob at each iteration
        g=gg(end);
        %     gg' % *** DEBUG ***
        m=m1;       % back up to previous iteration
        v=v1;
        w=w1;
        mm=sum(m,1)/k;
        f=(m(:)'*m(:)-k*mm(:)'*mm(:))/sum(v(:));
    end
    if ~fv
        m=m.*sx0(ones(k,1),:)+mx0(ones(k,1),:);	% unscale means
        v=v.*repmat(sx0.^2,k,1);                % and variances
    else
        v1=v;
        v=zeros(p,p,k);
        mk=eye(p)==1;                           % mask for diagonal elements
        v(repmat(mk,[1 1 k]))=v1';              % set from v1
    end
end
if fv              % check if full covariance matrices were requested
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Full Covariance matrices  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pl=p*(p+1)/2;
    lix=1:p^2;
    cix=repmat(1:p,p,1);
    rix=cix';
    lix(cix>rix)=[];                                        % index of lower triangular elements
    cix=cix(lix);                                           % index of lower triangular columns
    rix=rix(lix);                                           % index of lower triangular rows
    dix=find(rix==cix);
    lixi=zeros(p,p);
    lixi(lix)=1:pl;
    lixi=lixi';
    lixi(lix)=1:pl;                                        % reverse index to build full matrices
    v=reshape(v,p^2,k);
    v=v(lix,:)';                                            % lower triangular in rows

    % If data size is large then do calculations in chunks

    nb=min(n,max(1,floor(memsize/(24*p*k))));    % chunk size for testing data points
    nl=ceil(n/nb);                  % number of chunks
    jx0=n-(nl-1)*nb;                % size of first chunk
    %
    th=(l-floor(l))*n;
    sd=(nargout > 3*(l~=0)); % = 1 if we are outputting log likelihood values
    lp=floor(l)+sd;   % extra loop needed to calculate final G value
    %
    lpx=zeros(1,n);             % log probability of each data point
    wk=ones(k,1);
    wp=ones(1,p);
    wpl=ones(1,pl);             % 1 index for lower triangular matrix
    wnb=ones(1,nb);
    wnj=ones(1,jx0);

    % EM loop

    g=0;                        % dummy initial value for comparison
    gg=zeros(lp+1,1);
    ss=sd;                      % initialize stopping count (0 or 1)
    vi=zeros(p*k,p);            % stack of k inverse cov matrices each size p*p
    vim=zeros(p*k,1);       	% stack of k vectors of the form inv(v)*m
    mtk=vim;                  	% stack of k vectors of the form m
    lvm=zeros(k,1);
    wpk=repmat((1:p)',k,1);
    for j=1:lp
        g1=g;               	% save previous log likelihood (2*pi factor omitted)
        m1=m;                	% save previous means, variances and weights
        v1=v;
        w1=w;
        for ik=1:k

            % these lines added for debugging only
            %             vk=reshape(v(k,lixi),p,p);
            %             condk(ik)=cond(vk);
            %%%%%%%%%%%%%%%%%%%%
            [uvk,dvk]=eig(reshape(v(ik,lixi),p,p));	% convert lower triangular to full and find eigenvalues
            dvk=max(diag(dvk),c);                	% apply variance floor to eigenvalues
            vik=-0.5*uvk*diag(dvk.^(-1))*uvk';      % calculate inverse
            vi((ik-1)*p+(1:p),:)=vik;               % vi contains all mixture inverses stacked on top of each other
            vim((ik-1)*p+(1:p))=vik*m(ik,:)';       % vim contains vi*m for all mixtures stacked on top of each other
            mtk((ik-1)*p+(1:p))=m(ik,:)';           % mtk contains all mixture means stacked on top of each other
            lvm(ik)=log(w(ik))-0.5*sum(log(dvk));       % vm contains the weighted sqrt of det(vi) for each mixture
        end
        %
        %         % first do partial chunk
        %
        jx=jx0;
        ii=1:jx;
        xii=xs(ii,:).';
        py=reshape(sum(reshape((vi*xii-vim(:,wnj)).*(xii(wpk,:)-mtk(:,wnj)),p,jx*k),1),k,jx)+lvm(:,wnj);
        mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
        px=exp(py-mx(wk,:));            % find normalized probability of each mixture for each datapoint
        ps=sum(px,1);                   % total normalized likelihood of each data point
        px=px./ps(wk,:);                % relative mixture probabilities for each data point (columns sum to 1)
        lpx(ii)=log(ps)+mx;
        pk=px*wx(ii);                       % effective number of data points for each mixture (could be zero due to underflow)
        sx=px*xsw(ii,:);
        sx2=px*(xsw(ii,rix).*xs(ii,cix));	% accumulator for variance calculation (lower tri cov matrix as a row)
        for il=2:nl
            ix=jx+1;
            jx=jx+nb;        % increment upper limit
            ii=ix:jx;
            xii=xs(ii,:).';
            py=reshape(sum(reshape((vi*xii-vim(:,wnb)).*(xii(wpk,:)-mtk(:,wnb)),p,nb*k),1),k,nb)+lvm(:,wnb);
            mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
            px=exp(py-mx(wk,:));            % find normalized probability of each mixture for each datapoint
            ps=sum(px,1);                   % total normalized likelihood of each data point
            px=px./ps(wk,:);                % relative mixture probabilities for each data point (columns sum to 1)
            lpx(ii)=log(ps)+mx;
            pk=pk+px*wx(ii);                    % effective number of data points for each mixture (could be zero due to underflow)
            sx=sx+px*xsw(ii,:);             % accumulator for mean calculation
            sx2=sx2+px*(xsw(ii,rix).*xs(ii,cix));	% accumulator for variance calculation
        end
        g=lpx*wx;                   % total log probability summed over all data points
        gg(j)=g;                    % save convergence history
        w=pk/nwt;               	% w(k,1) normalize to get the column of weights
        if pk                       % if all elements of pk are non-zero
            m=sx./pk(:,wp);         % find mean and mean square
            v=sx2./pk(:,wpl);
        else
            wm=pk==0;                       % mask indicating mixtures with zero weights
            nz=sum(wm);                  % number of zero-weight mixtures
            [vv,mk]=sort(lpx);             % find the lowest probability data points
            m=zeros(k,p);                   % initialize means and variances to zero (variances are floored later)
            v=zeros(k,pl);
            m(wm,:)=xs(mk(1:nz),:);                % set zero-weight mixture means to worst-fitted data points
            w(wm)=1/n;                      % set these weights non-zero
            w=w*n/(n+nz);                   % normalize so the weights sum to unity
            wm=~wm;                         % mask for non-zero weights
            m(wm,:)=sx(wm,:)./pk(wm,wp);  % recalculate means and variances for mixtures with a non-zero weight
            v(wm,:)=sx2(wm,:)./pk(wm,wpl);
        end
        v=v-m(:,cix).*m(:,rix);                 % subtract off mean squared
        if g-g1<=th && j>1
            if ~ss, break; end  %  stop
            ss=ss-1;       % stop next time
        end
    end
    if sd  % we need to calculate the final probabilities
        pp=lpx'-0.5*p*log(2*pi)-lsx;   % log of total probability of each data point
        gg=gg(1:j)/nwt-0.5*p*log(2*pi)-lsx;    % average log prob at each iteration
        g=gg(end);
        %             gg' % *** DEBUG ONLY ***
        m=m1;                                           % back up to previous iteration
        v=zeros(p,p,k);                                 % reserve spave for k full covariance matrices
        trv=0;                                          % sum of variance matrix traces
        for ik=1:k                                      % loop for each mixture to apply variance floor
            [uvk,dvk]=eig(reshape(v1(ik,lixi),p,p));	% convert lower triangular to full and find eigenvectors
            dvk=max(diag(dvk),c);                       % apply variance floor to eigenvalues
            v(:,:,ik)=uvk*diag(dvk)*uvk';               % reconstitute full matrix
            trv=trv+sum(dvk);                           % add trace to the sum
        end
        w=w1;
        mm=sum(m,1)/k;
        f=(m(:)'*m(:)-k*mm(:)'*mm(:))/trv;
    else
        v1=v;                                           % lower triangular form
        v=zeros(p,p,k);                                 % reserve spave for k full covariance matrices
        for ik=1:k                                      % loop for each mixture to apply variance floor
            [uvk,dvk,]=eig(reshape(v1(ik,lixi),p,p));	% convert lower triangular to full and find eigenvectors
            dvk=max(diag(dvk),c);                       % apply variance floor
            v(:,:,ik)=uvk*diag(dvk)*uvk';               % reconstitute full matrix
        end
    end
    m=m.*sx0(ones(k,1),:)+mx0(ones(k,1),:);  % unscale means
    v=v.*repmat(sx0'*sx0,[1 1 k]);
end
if l==0         % suppress the first three output arguments if l==0
    m=g;
    v=f;
    w=pp;
end

