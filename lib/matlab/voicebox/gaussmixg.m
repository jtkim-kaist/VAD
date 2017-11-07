function [mg,vg,pg,pv]=gaussmixg(m,v,w,n)
%GAUSSMIXG global mean, variance and mode of a GMM
%
% Usage: (1) gaussmixg(m,v,w)               % plot the mean and mode positions of a GMM
%        (2) [mg,vg]=gaussmixg(m,v,w)       % find global mean and covariance of a GMM
%        (3) [mg,vg,pg]=gaussmixg(m,v,w)    % find global mean,covariance and mode of a GMM
%        (4) [mg,vg,pg,pv]=gaussmixg(m,v,w) % ... also find log probability of the peak
%
%  Inputs:  M(k,p) = mixture means for pg(p)
%           V(k,p) or V(p,p,k) variances (diagonal or full)
%           W(k,1) = mixture weights
%           N      = maximum number of modes to find [default 1]
%
% Outputs: MG(1,p) = global mean
%          VG(p,p) = global covariance
%          PG(N,p) = sorted list of N modes
%          PV(N,1) = log pdf at the modes PG(N,p) (in decreasing order)
%
%  This routine finds the global mean and covariance matrix of a Gaussian Mixture (GMM). It also
%  attempts to find up to N local maxima using a combination of the fixed point and quadratic
%  Newton-Raphson algorithms from [1]. Currently, N must be less than or equal to the number of
%  mixtures K. In general the PDF surface of a GMM can be very complicated with many local maxima [2]
%  and, as discussed in [1,2], this algorithm is not guaranteed to find the N highest. In [2], it is
%  conjectured that the number of local maxima is <=K for the following cases (a) P=1, (b) all
%  mixture covariance matrices are equal and (c) all mixture covariance matrices are multiples of
%  the identity.
%
% Refs:
%   [1]	M. Á. Carreira-Perpiñán. Mode-finding for mixtures of gaussian distributions.
%       IEEE Trans. Pattern Anal and Machine Intell, 22 (11): 1318–1323, 2000. doi: 10.1109/34.888716.
%   [2] M. Á. Carreira-Perpiñán and C. K. I. Williams. On the number of modes of a gaussian mixture.
%       In Proc Intl Conf on Scale Space Theories in Computer Vision, volume LNCS 2695, pages 625–640,
%       Isle of Skye, June 2003. doi: 10.1007/3-540-44935-3_44.

% Bugs/Suggestions:
% (1) Sometimes the mode is not found, e.g. m=[0 1; 1 0];v=[.01 10; 10 .01];
%     has a true mode near (0,0). Could add to the list of mode candidates
%     all the pairwise intersections of the mixtures.
%     Another is: m=[0 0; 10 0.3]; v=[1 1; 1000 .001];
% (2) When merging candidates, we should keep the one with the highest probability
% (3) could preserve the fixed arrays between calls if p and/or k are unchanged
%
% See also: gaussmix, gaussmixd, gaussmixp, randvec

%      Copyright (C) Mike Brookes 2000-2012
%      Version: $Id: gaussmixg.m 3227 2013-07-04 15:42:04Z dmb $
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
% Algorithm parameters
nfp=2;          % number of fixed point iterations to do at start
maxloop=60;     % maximum number of iterations
ssf=0.1;        % factor when calculating minimum mode separation
%Sort out arguments
[k,p]=size(m);
if nargin<4
    n=1;
    if nargin<3
        w=ones(k,1);
        if nargin<2
            v=ones(k,p);
        end
    end
end
if ~nargout
    if nargin<4
        n=k;
    end
    nao=4;
else
    nao=nargout;
end

full=numel(size(v))>2 || k==1 && numel(v)>p; % test for full covariance matrices
if full && p==1   % if p=1 then we force diagonal covariance matrices
    v=reshape(v,1,k)';
    full=false;
end
w=w/sum(w);  % force w to sum to unity
% calculate the global mean and covariance
mg=w.'*m;
mz=m-mg(ones(k,1),:);   % means relative to global mean
if nao>2                % need to calculate the modes
    nx=k;               % number of pg values initially
    kk=reshape(repmat(1:nx,k,1),k*nx,1); % [1 1 1 2 2 2 ... nx nx nx]'
    km=repmat(1:k,1,nx)'; % [1:k 1:k ... 1:k]'
    % sort out indexing for all data value pairs; needed to eliminate duplicates
    nxp=nx*(nx-1)/2;
    ja=(1:nxp)';
    jb=floor((3+sqrt(8*ja-3))/2);       % [2 3 3 4 4 4 ... nx]'
    ja=ja-(jb.^2-3*jb+2)/2;             % [1 1:2 1:3 ... 1:nx-1]'
    jc=ones(nxp,1);
    % sort out indexing for vectorized upper triagular matrix
    npu=p*(p+1)/2;                      % number of distinct elements in a symmetrial (p,p) matrix
    kw=1:npu;
    ku=floor((1+sqrt(8*kw-3))/2);       % [1 2 2 3 3 3 ... p]
    kv=kw-(ku.^2-ku)/2;                 % [1 1:2 1:3 ... 1:p]
    zpp=zeros(p,p);
    zpp(kv+p*(ku-1))=kw;
    zpp(ku+p*(kv-1))=kw;
    kw=reshape(zpp,1,[]);               % maps vectorized upper triangular to vectorized full matrix
    kp=repmat(1:p,1,p);                 % row indices for a (p,p) matrix as a row vector
    kq=reshape(repmat(1:p,p,1),1,p^2);  % col indices for a (p,p) matrix as a row vector
    kr=p*kp-p+kq;                       % transpose indexing for a vectorized (p,p) matrix
    kd=1:p+1:p^2;                       % diagonal indices of a (p,p) matrix
    % unity vectors to make efficient replication
    wk=ones(k,1);
    wnx=ones(nx,1);
    if full
        vg=mz.'*(mz.*w(:,ones(1,p)))+reshape(reshape(v,p^2,k)*w,p,p);
        % now determine the mode
        vi=zeros(p*k,p);                    % stack of k inverse cov matrices each size p*p times -0.5
        vim=zeros(p*k,1);                   % stack of k vectors of the form -0.5*inv(vt)*m
        mtk=vim;                             % stack of k vectors of the form m
        lvm=zeros(k,1);
        wpk=repmat((1:p)',k,1);
        for i=1:k    % loop for each mixture
            [uvk,dvk]=eig(v(:,:,i));      % find eigenvalues
            dvk=diag(dvk);
            if(any(dvk<=0))
                error('Covariance matrix for mixture %d is not positive definite',i);
            end
            vik=-0.5*uvk*diag(dvk.^(-1))*uvk';   % calculate inverse including -0.5 factor
            vi((i-1)*p+(1:p),:)=vik;           % vi contains all mixture inverses stacked on top of each other
            vim((i-1)*p+(1:p))=vik*m(i,:)';   % vim contains vi*m for all mixtures stacked on top of each other
            mtk((i-1)*p+(1:p))=m(i,:)';       % mtk contains all mixture means stacked on top of each other
            lvm(i)=log(w(i))-0.5*sum(log(dvk));       % vm contains the weighted sqrt of det(vi) for each mixture
        end
        vif=reshape(permute(reshape(vi,p,k,p),[2 1 3]),k,p^2); % each covariance matrix as a vectorized row
        vimf=reshape(vim,p,k)'; % vi*m as a row for each mixture
        ss=sqrt(min(v(repmat(kd,k,1)+repmat(p^2*(0:k-1)',1,p)),[],1))*ssf/sqrt(p);   % minimum separation of modes [this is a conservative guess]
    else
        vg=mz.'*(mz.*w(:,ones(1,p)))+diag(w.'*v);
        % now determine the mode
        vi=-0.5*v.^(-1);                % vi(k,p) = data-independent scale factor in exponent
        vi2=vi(:,ku).*vi(:,kv);         % vi2(k,npu) = upper triangular Hessian data dependent term
        lvm=log(w)-0.5*sum(log(v),2);   % log of external scale factor (excluding -0.5*p*log(2pi) term)
        vim=vi.*m;
        vim2=vim(:,ku).*vim(:,kv);      % vim2(k,npu) = upper triangular Hessian data independent term
        vimvi=vim(:,kp).*vi(:,kq);      % vimvi(k,p^2) = vectorized Hessian term
        ss=sqrt(min(v,[],1))*ssf/sqrt(p);   % minimum separation of modes [this is a conservative guess]
    end
    pgf=zeros(nx,p);                    % space for fixed point update
    sv=0.01*ss;                         % convergence threshold
    pg=m;  % initialize mode candidates to mixture means
    i=1;   % loop counter
    %     gx=zeros(nx,p);   %%%%%%%%%%%% temp
    while i<=maxloop
        %         pg00=pg0;   %%%%%%%%%%%% temp
        pg0=pg;                         % save previous mode candidates pg(nx,p)
        if full
            py=reshape(sum(reshape((vi*pg'-vim(:,wnx)).*(pg(:,wpk)'-mtk(:,wnx)),p,nx*k),1),k,nx)+lvm(:,wnx);
            mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
            px=exp(py-mx(wk,:));            % find normalized probability of each mixture for each datapoint
            ps=sum(px,1);                   % total normalized likelihood of each data point
            px=(px./ps(wk,:))';             % px(nx,k) = relative mixture probabilities for each data point (rows sum to 1)
            % calculate the fixed point update
            pxvif=px*vif;    % pxvif(nx,p^2)
            pxvimf=px*vimf;  % pxvimf(nx,p)
            for j=1:nx
                pgf(j,:)=pxvimf(j,:)/reshape(pxvif(j,:),p,p);
            end
        else
            py=reshape(sum((pg(kk,:)-m(km,:)).^2.*vi(km,:),2),k,nx)+lvm(:,wnx);
            mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
            px=exp(py-mx(wk,:));            % find normalized probability of each mixture for each datapoint
            ps=sum(px,1);                   % total normalized likelihood of each data point
            px=(px./ps(wk,:))';             % px(nx,k) = relative mixture probabilities for each data point (rows sum to 1)
            % calculate the fixed point update
            pxvim=px*vim;
            pxvi=px*vi;
            pgf=pxvim./pxvi;       % fixed point update for all points
        end
        if i>nfp
            % calculate gradient and Hessian; see [1] equations (4) and (5)
            lp=log(ps)+mx;              % log prob of each data point
            if full
                %                 gx0=gx;   %%%%%%%%%%%% temp
                gx=pxvimf-reshape(sum(repmat(pg,p,1).*reshape(pxvif,[],p),2),nx,p);
                vimpg=repmat(vimf,nx,1)-reshape(permute(reshape(pg*vi',nx,p,k),[3 1 2]),[],p); % vimpg(k*nx,p)
                hx1=2*reshape(sum(reshape(repmat(reshape(px',[],1),1,npu).*vimpg(:,ku).*vimpg(:,kv),k,[]),1),nx,[]);
                hx=pxvif+hx1(:,kw);
            else
                gx=pxvim-pxvi.*pg;               % gradient for each data point (one row per point)
                hx1=px*vim2+(px*vi2).*pg(:,ku).*pg(:,kv);
                hx2=(px*(vimvi)).*pg(:,kq);
                hx=2*(hx1(:,kw)-hx2-hx2(:,kr));
                hx(:,kd)=hx(:,kd)+pxvi;
            end
            hx=reshape(hx',p,p,nx);
            for j=1:nx
                if all(eig(hx(:,:,j))<0)  % if positive definite
                    pg(j,:)=pg(j,:)+gx(j,:)/hx(:,:,j); % do a Newton-Raphson update
                    if full
                        pyj=sum(reshape((vi*pg(j,:)'-vim).*(pg(j,wpk)'-mtk),p,k),1)'+lvm;
                    else
                        pyj=sum((repmat(pg(j,:),k,1)-m).^2.*vi,2)+lvm;
                    end
                    mxj=max(pyj);                % find normalizing factor for each data point to prevent underflow when using exp()
                    pxj=exp(pyj-mxj);            % find normalized probability of each mixture for each datapoint
                    psj=sum(pxj,1);                   % total normalized likelihood of each data point
                    lpj=log(psj)+mxj;              % log prob of updated data point
                    if lpj<lp(j)       % check if the probability has decreased
                        pg(j,:)=pgf(j,:);   % if so, do fixed point update
                    end
                else
                    pg(j,:)=pgf(j,:);   % else do fixed point update
                end
            end
        else
            pg=pgf;       % fixed point update for all points
        end
        if all(all(abs(pg-pg0)<sv(wnx,:))) && i+2<maxloop
            maxloop=min(maxloop,i+2);        % two more loops if converged if converged
        end
        %         [all(pg==pgf,2) pg [i; repmat(NaN,nx-1,1)] pg-pg0]  %debug: [fixed-point x-mode iteration delta-x]
        jd=all(abs(pg(jb,:)-pg(ja,:))<ss(jc,:),2);   % find duplicate modes
        if any(jd)
            jx=sparse([(1:nx)';ja;jb],[(1:nx)';jb;ja],[wnx;jd;jd]);   % neighbour matrix
            kx=any((jx*jx)>0 & ~jx,2);  % find chains that  are not fully connected
            while any(kx)
                kx=any(jx(:,kx),2);     % destroy all links connected to these chains
                jx(kx,:)=0;
                jx(:,kx)=0;
                % jx(kx,kx)=1;
                kx=any((jx*jx)>0 & ~jx,2);
            end
            jx([1:nx+1:nx*nx (ja+(jb-1)*nx)'])=0; % reset the upper triangle + diagonal
            pg(any(jx,2),:)=[];   % delete the duplicates
            % update nx and anything that depends on it
            nx=size(pg,1);
            pgf=zeros(nx,p);                    % space for fixed point update
            wnx=ones(nx,1);
            nxp=nx*(nx-1)/2;
            ja=ja(1:nxp);
            jb=jb(1:nxp);
            kk=reshape(repmat(1:nx,k,1),k*nx,1); % [1 1 1 2 2 2 ... nx nx nx]'
            km=reshape(repmat(1:k,1,nx),k*nx,1); % [1:k 1:k ... 1:k]'
            jc=ones(nxp,1);
        end
        i=i+1;
    end
    %     calculate the log pdf at each mode
    if full
        py=reshape(sum(reshape((vi*pg'-vim(:,wnx)).*(pg(:,wpk)'-mtk(:,wnx)),p,nx*k),1),k,nx)+lvm(:,wnx);
        mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
        px=exp(py-mx(wk,:));            % find normalized probability of each mixture for each datapoint
        ps=sum(px,1);                   % total normalized likelihood of each data point
    else
        py=reshape(sum((pg(kk,:)-m(km,:)).^2.*vi(km,:),2),k,nx)+lvm(:,wnx);
        mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
        px=exp(py-mx(wk,:));            % find normalized probability of each mixture for each datapoint
        ps=sum(px,1);                   % total normalized likelihood of each data point
    end
    [pv,ix]=sort((log(ps)+mx)'-0.5*p*log(2*pi),'descend');
    pg=pg(ix,:);
    if n<numel(pv) % only keep the first n modes
        pg=pg(1:n,:);
        pv=pv(1:n);
    end
elseif nao>1
    if full
        vg=mz.'*(mz.*w(:,ones(1,p)))+reshape(reshape(v,p^2,k)*w,p,p);
    else
        vg=mz.'*(mz.*w(:,ones(1,p)))+diag(w.'*v);
    end
end

if ~nargout
    % now plot the result
    clf;
    pg1=pg(1,:);
    lpm=gaussmixp(mg,m,v,w);
    lpp=pv(1);
    switch p
        case 1
            gaussmixp([],m,v,w);
            hold on;
            ylim=get(gca,'ylim')';
            plot([mg mg]',ylim,'-k',[mg mg; mg mg]+[-1 1; -1 1]*sqrt(vg),ylim,':k');
            plot(pg(1),pv(1),'^k');
            if numel(pg)>1
                plot(pg(2:end),pv(2:end),'xk');
            end
            hold off;
            title(sprintf('Mean+-sd = %.3g+-%.3g LogP = %.3g, Mode\\Delta = %.3g LogP = %.3g',mg,sqrt(vg),lpm,pg1,lpp));
            xlabel('x');
        case 2
            gaussmixp([],m,v,w);
            hold on;
            t=linspace(0,2*pi,100);
            xysd=chol(vg)'*[cos(t); sin(t)]+repmat(mg',1,length(t));
            plot(xysd(1,:),xysd(2,:),':k',mg(1),mg(2),'ok');
            plot(pg(1,1),pg(1,2),'^k');
            if numel(pv)>1
                plot(pg(2:end,1),pg(2:end,2),'xk');
            end
            hold off;
            title(sprintf('Mean = (%.3g,%.3g) LogP = %.3g, Mode:\\Delta = (%.3g,%.3g) LogP = %.3g',mg,lpm,pg1,pv(1)));
            xlabel('x');
            ylabel('y');
        otherwise
            nx=200;
            nc=ceil(sqrt(p/2));
            nr=ceil(p/nc);
            sdx=sqrt(diag(vg))';  % std deviation
            minx=min([mg; pg],[],1)-1.5*sdx;
            maxx=max([mg; pg],[],1)+1.5*sdx;
            ix=2:p; % selected indices
            for i=1:p
                xi=linspace(minx(i),maxx(i),nx)';
                [mm,vm,wm]=gaussmixd(mg(ix),m,v,w,[],ix);
                ym=gaussmixp(xi,mm,vm,wm')+lpm-gaussmixp(mg(i),mm,vm,wm');
                [mp,vp,wp]=gaussmixd(pg(1,ix),m,v,w,[],ix);
                yp=gaussmixp(xi,mp,vp,wp')+lpp-gaussmixp(pg1(i),mp,vp,wp');
                subplot(nr,nc,i);
                plot(xi,ym,'-k',mg(i),lpm,'ok',xi,yp,':b',pg1(i),lpp,'^b');
                axisenlarge([-1 -1 -1 -1.05]);
                hold on
                plot([mg(i) mg(i)],get(gca,'ylim'),'-k');
                hold off
                xlabel(sprintf('x[%d], Mean = %.3g, Mode\\Delta = %.3g',i,mg(i),pg1(i)));
                if i==1
                    title(sprintf('Log Prob: Mean = %.3g, Mode\\Delta = %.3g',lpm,lpp));
                end
                ix(i)=i;
            end
    end
end
