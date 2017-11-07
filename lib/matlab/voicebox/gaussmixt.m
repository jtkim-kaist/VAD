function [m,v,w]=gaussmixt(m1,v1,w1,m2,v2,w2)
%GAUSSMIXT Multiply two GMM pdfs
%
% Inputs: Input mixtures: k1,k2 mixtures, p dimensions
%
%   M(k1,p) = mixture means for mixture 1
%   V(k1,p) or V(p,p,k1) variances (diagonal or full)
%   W(k1,1) = mixture weights
%   M(k2,p) = mixture means for mixture 2
%   V(k2,p) or V(p,p,k2) variances (diagonal or full)
%   W(k2,1) = mixture weights
%
% Outputs:
%
%   M(k1*k2,p) = mixture means
%   V(k1*k2,p) or V(p,p,k1*k2) if p>1 and at least one input has full covariance matrix
%   W(k1*k2,1) = mixture weights
%
% See also: gaussmix, gaussmixg, gaussmixp, randvec

%      Copyright (C) Mike Brookes 2000-2012
%      Version: $Id: gaussmixt.m 5453 2014-11-19 13:10:51Z dmb $
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
persistent r13 r21 r22 r31 r312 r112 r1223 r321 ch1h r122 r124
if isempty(r21)
    r13=[1 3];
    r21=[2 1];
    r22=[2 2];
    r31=[3 1];
    r112=[1 1 2];
    r122=[1 2 2];
    r124=[1 2 4];
    r312=[3 1 2];
    r321=[3 2 1];
    r1223=[1 2 2 3];
    ch1h=[-0.5; 1; -0.5];
end
[k1,p]=size(m1);
[k2,p2]=size(m2);
f1=ndims(v1)>2 || size(v1,1)>k1; % full covariance matrix is supplied
f2=ndims(v2)>2 || size(v2,1)>k2; % full covariance matrix is supplied
% ff=f1+2*f2;
if p~=p2
    error('mixtures must have the same vector dimension');
end
k=k1*k2;
j1=repmat((1:k1)',k2,1);
j2=reshape(repmat(1:k2,k1,1),k,1);
if p==1
    % display('1D vectors');
    p1=1./v1(:);
    p2=1./v2(:);
    v=1./(p1(j1)+p2(j2));
    s1=p1.*m1;
    s2=p2.*m2;
    m=v.*(s1(j1)+s2(j2));
    v12=v1(j1)+v2(j2);
    wx=-0.5*(m1(j1)-m2(j2)).^2./v12(:);
    wx=wx-max(wx); % normalize to avoid underflow
    w=w1(j1).*w2(j2).*exp(wx)./sqrt(v12(:));
    w=w/sum(w);
else
    if ~f1 && ~f2 % both diagonal covariances
        % display('both diagonal');
        p1=1./v1;
        p2=1./v2;
        v=1./(p1(j1,:)+p2(j2,:));
        s1=p1.*m1;
        s2=p2.*m2;
        m=v.*(s1(j1,:)+s2(j2,:));
        v12=v1(j1,:)+v2(j2,:);
        wx=-0.5*sum((m1(j1,:)-m2(j2,:)).^2./v12,2);
        wx=wx-max(wx); % normalize to avoid underflow
        w=w1(j1).*w2(j2).*exp(wx)./sqrt(prod(v12,2));
        w=w/sum(w);
    else % at least one full covariances
        m=zeros(k,p);
        v=zeros(p,p,k);
        w=zeros(k,1);
        wx=w;
        idp=1:p+1:p*p; % diagonal elements of p x p matrix
        if p==2                 % special code for 2D vectors
            if ~f2  % GMM 2 is diagonal
                % display('2D GMM 2 diagonal');
                p2=1./v2;
                pm2=p2.*m2;
                vx1=permute(v1,r312);
                vx1=vx1(:,r124);
                px1=vx1./repmat((vx1(:,1).*vx1(:,3)-vx1(:,2).^2),1,3); % [a b; b c] -> [c -b a]
                pm1=m1.*px1(:,r31)-m1(:,r21).*px1(:,r22);
                px=px1(j1,:);
                px(:,r31)=px(:,r31)+p2(j2,:);  % add onto diagonal elements
                vijx=vx1(j1,:);
                vijx(:,r13)=vijx(:,r13)+v2(j2,:);  % add onto diagonal elements
            elseif ~f1 % GMM 1 is diagonal
                % display('2D GMM 1 diagonal');
                p1=1./v1;
                pm1=p1.*m1;
                vx2=permute(v2,r312);
                vx2=vx2(:,r124);
                px2=vx2./repmat((vx2(:,1).*vx2(:,3)-vx2(:,2).^2),1,3); % [a b; b c] -> [c -b a]
                pm2=m2.*px2(:,r31)-m2(:,r21).*px2(:,r22);
                px=px2(j2,:);
                px(:,r31)=px(:,r31)+p1(j1,:);  % add onto diagonal elements
                vijx=vx2(j2,:);
                vijx(:,r13)=vijx(:,r13)+v1(j1,:);  % add onto diagonal elements
            else % both full covariances
                % display('2D both full');
                vx1=permute(v1,r312);
                vx1=vx1(:,r124); % make each 2 x 2 matrix into a row [a b; b c] -> [a b c]
                px1=vx1./repmat((vx1(:,1).*vx1(:,3)-vx1(:,2).^2),1,3); % [a b; b c] -> [c -b a]
                vx2=permute(v2,r312);
                vx2=vx2(:,r124);
                px2=vx2./repmat((vx2(:,1).*vx2(:,3)-vx2(:,2).^2),1,3); % [a b; b c] -> [c -b a]
                pm1=m1.*px1(:,r31)-m1(:,r21).*px1(:,r22);
                pm2=m2.*px2(:,r31)-m2(:,r21).*px2(:,r22);
                px=px1(j1,:)+px2(j2,:);
                vijx=vx1(j1,:)+vx2(j2,:);
            end
            vx=px./repmat((px(:,1).*px(:,3)-px(:,2).^2),1,3);   % divide by determinant to get inverse
            m=pm1(j1,:)+pm2(j2,:);
            m=m.*vx(:,r13)+m(:,r21).*vx(:,r22);                 % multiple by 2 x 2 matrix vx
            v=reshape(vx(:,r1223)',[2 2 k]);                    % convert vx to a 3D array of 2 x 2 matrices
            m12=m1(j1,:)-m2(j2,:);                              % subtract means to calculate weight exponent
            dij=vijx(:,1).*vijx(:,3)-vijx(:,2).^2;              % determinant of V1+V2
            wx=m12(:,r112).*m12(:,r122).*vijx(:,r321)*ch1h./dij;% exponent of weight
            w=w1(j1).*w2(j2)./sqrt(dij);                        % weight is w*exp(wx)
        else
            if ~f2  % GMM 2 is diagonal
                % display('GMM 2 diagonal');
                p2=1./v2;
                pm2=p2.*m2;
                for i=1:k1
                    v1i=v1(:,:,i);
                    p1i=inv(v1i);
                    m1i=m1(i,:);
                    pm1i=m1i*p1i;
                    w1i=w1(i);
                    ix=i;
                    for j=1:k2
                        pij=p1i;
                        pij(idp)=pij(idp)+p2(j,:);
                        vix=inv(pij);
                        vij=v1i;
                        vij(idp)=vij(idp)+v2(j,:);
                        v(:,:,ix)=vix;
                        m(ix,:)=(pm2(j,:)+pm1i)*vix;
                        m12=m2(j,:)-m1i;
                        wx(ix)=-0.5*m12/vij*m12';           % exponent of weight
                        w(ix)=w2(j)*w1i/sqrt(det(vij));     % weight is w*exp(wx)
                        ix=ix+k1;
                    end
                end
            elseif ~f1 % GMM 1 is diagonal
                % display('GMM 1 diagonal');
                p1=1./v1;
                pm1=p1.*m1;
                ix=1;
                for j=1:k2
                    v2j=v2(:,:,j);
                    p2j=inv(v2j);
                    m2j=m2(j,:);
                    pm2j=m2j*p2j;
                    w2j=w2(j);
                    for i=1:k1
                        pij=p2j;
                        pij(idp)=pij(idp)+p1(i,:);
                        vix=inv(pij);
                        vij=v2j;
                        vij(idp)=vij(idp)+v1(i,:);
                        v(:,:,ix)=vix;
                        m(ix,:)=(pm1(i,:)+pm2j)*vix;
                        m12=m1(i,:)-m2j;
                        wx(ix)=-0.5*m12/vij*m12';           % exponent of weight
                        w(ix)=w1(i)*w2j/sqrt(det(vij));     % weight is w*exp(wx)
                        ix=ix+1;
                    end
                end
            else % both full covariances
                % display('both full');
                p1=zeros(p,p,k1);
                pm1=zeros(k1,p);
                for i=1:k1
                    p1i=inv(v1(:,:,i));
                    p1(:,:,i)=p1i;
                    pm1(i,:)=m1(i,:)*p1i;
                end
                ix=1;
                for j=1:k2
                    v2j=v2(:,:,j);
                    p2j=inv(v2j);
                    m2j=m2(j,:);
                    pm2j=m2j*p2j;
                    w2j=w2(j);
                    for i=1:k1
                        pij=p1(:,:,i)+p2j;
                        vix=inv(pij);
                        v(:,:,ix)=vix;
                        vij=v1(:,:,i)+v2j;
                        m(ix,:)=(pm1(i,:)+pm2j)*vix;
                        m12=m1(i,:)-m2j;
                        wx(ix)=-0.5*m12/vij*m12';           % exponent of weight
                        w(ix)=w1(i)*w2j/sqrt(det(vij));     % weight is w*exp(wx)
                        ix=ix+1;
                    end
                end
            end
            
        end
        wx=wx-max(wx);              % adjust exponents to avoid underflow
        w=w.*exp(wx);               % calculate weights
        w=w/sum(w);                 % normalize weights to sum to unity
        if k==1
            v=reshape(v,size(v,1),size(v,2)); % squeeze last dimension of v if possible
        end
    end
end
if ~nargout
    if p==1
        nxx=256; % number of points to plot
        nsd=3; % number of std deviations
        sd=sqrt([v1(:);v2(:);v]);
        ma=[m1;m2;m];
        xax=linspace(min(ma-nsd*sd),max(ma+nsd*sd),nxx);
        plot(xax,gaussmixp(xax(:),m1,v1,w1),'--b');
        hold on
        plot(xax,gaussmixp(xax(:),m2,v2,w2),':r');
        plot(xax,gaussmixp(xax(:),m,v,w),'-k');
        hold off
        ylabel('Log probability density');
        legend('Mix 1','Mix 2','Product','location','best');
        axisenlarge([-1 -1 -1 -1.05]);
    elseif p==2
        nxx=128; % number of points to plot
        nsd=3;
        if f1
            s1=sqrt([v1(1:4:end)' v1(4:4:end)']); % extract diagonal elements only
        else
            s1=sqrt(v1);
        end
        if f2
            s2=sqrt([v2(1:4:end)' v2(4:4:end)']); % extract diagonal elements only
        else
            s2=sqrt(v2);
        end
        if ndims(v)>2 || size(v,1)>k
            s3=sqrt([v(1:4:end)' v(4:4:end)']); % extract diagonal elements only
        else
            s3=sqrt(v);
        end
        mal=[m1;m2;m];
        sal=[s1;s2;s3];
        xax=linspace(min(mal(:,1)-nsd*sal(:,1)),max(mal(:,1)+nsd*sal(:,1)),nxx);
        yax=linspace(min(mal(:,2)-nsd*sal(:,2)),max(mal(:,2)+nsd*sal(:,2)),nxx);
        xx(:,:,1)=repmat(xax',1,nxx);
        xx(:,:,2)=repmat(yax,nxx,1);
        xx=reshape(xx,nxx^2,2);
        subplot(2,2,1);
        imagesc(xax,yax,reshape(gaussmixp(xx,m1,v1,w1),nxx,nxx)');
        axis 'xy';
        title('Input Mix 1');
        subplot(2,2,2);
        imagesc(xax,yax,reshape(gaussmixp(xx,m2,v2,w2),nxx,nxx)');
        axis 'xy';
        title('Input Mix 2');
        subplot(2,2,3);
        imagesc(xax,yax,reshape(gaussmixp(xx,m,v,w),nxx,nxx)');
        axis 'xy';
        title('Product GMM');
    end
end
