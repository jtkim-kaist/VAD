function [g,u,k]=gausprod(m,c,e)
%GAUSPROD calculates a product of gaussians [G,U,K]=(M,C)
% calculates the product of n d-dimensional multivariate gaussians
% this product is itself a gaussian
% Inputs: m(d,n) - each column is the mean of one of the gaussians
%         c(d,d,n) - contains the d#d covariance matrix for each gaussian
%                    Alternatives: (i) c(d,n) if diagonal (ii) c(n) if c*I or (iii) omitted if I
% not yet implemented: 
%         e(d,d,n) - contains orthogonal eigenvalue matrices and c(d,n) contains eigenvalues.
%                    Covariance matrix is E(:,:,k)*diag(c(:,k))*E(:,:,k)'
%                    c(d,n) can then contain 0 and Inf entries
%
% Outputs: g log gain (= log(integral) )
%          u(d,1) mean vector
%          k(d,d) or k(d) or k(1) = covariance matrix, diagonal or multiple of I (same form as input)
%

% Bugs: this version works with singular covariance matrices provided that their null spaces are distinct
%       we could improve it slightly by doing the pseudo inverses locally and keeping track of null spaces

%	   Copyright (C) Mike Brookes 2004
%      Version: $Id: gausprod.m 4777 2014-07-01 19:40:53Z dmb $
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

[d,n]=size(m);
if nargin>2
    error('third argument not yet implemented in gausprod');
end
if nargin<2     % all covariance matrices = I
    c=ones(n,1);
end
if ~nargout     % save input data for plotting
    m0=m;
    c0=c;
end

sc=size(c);
if length(sc)<3
    if(sc(2)==1) & (n>1)    % covariance matrices are multiples of the identity
        jj=1;
        jj2=2;
        gj=zeros(n,1);
        while jj<n
            for j=1+jj:jj2:n        % we combine the gaussians in pairs
                k=j-jj;
                cjk=c(k)+c(j);
                dm=m(:,k)-m(:,j);
                gj(k)=gj(k)+gj(j)-0.5*(d*log(cjk)+dm'*dm/cjk);
                m(:,k)=(c(k)*m(:,j)+c(j)*m(:,k))/cjk;
                c(k)=c(k)*c(j)/cjk;
            end
            jj=jj2;
            jj2=2*jj;
        end
        g=gj(1);
        k=c(1);
        u=m(:,1);
    else                    % diagonal covariance matrices
        jj=1;
        jj2=2;
        gj=zeros(n,1);
        while jj<n
            for j=1+jj:jj2:n        % we combine the gaussians in pairs
                k=j-jj;
                cjk=c(:,k)+c(:,j);
                dm=m(:,k)-m(:,j);
                ix=cjk>d*max(cjk)*eps;      % calculate the psedo inverse directly
                piv=zeros(d,1);
                piv(ix)=cjk(ix).^(-1);
                gj(k)=gj(k)+gj(j)-0.5*(log(prod(cjk))+piv'*dm.^2);
                m(:,k)=piv.*(c(:,k).*m(:,j)+c(:,j).*m(:,k));
                c(:,k)=c(:,k).*piv.*c(:,j);
            end
            jj=jj2;
            jj2=2*jj;
        end
        g=gj(1);
        k=c(:,1);
        u=m(:,1);
    end
else                        % full covariance matrices
    jj=1;
    jj2=2;
    gj=zeros(n,1);
    while jj<n
        for j=1+jj:jj2:n        % we combine the gaussians in pairs
            k=j-jj;
            cjk=c(:,:,k)+c(:,:,j);
            dm=m(:,k)-m(:,j);
            piv=pinv(cjk);
            gj(k)=gj(k)+gj(j)-0.5*(log(det(cjk))+dm'*piv*dm);
            m(:,k)=c(:,:,k)*piv*m(:,j)+c(:,:,j)*piv*m(:,k);
            c(:,:,k)=c(:,:,k)*piv*c(:,:,j);
            c(:,:,k)=0.5*(c(:,:,k)+c(:,:,k)');  % ensure exactly symmetric
        end
        jj=jj2;
        jj2=2*jj;
    end
    g=gj(1);
    k=c(:,:,1);
    u=m(:,1);
end
g=g-0.5*(n-1)*d*log(2*pi);

if ~nargout                 % plot results if no output arguments
    if d==1                 % one-dimensional vectors
        x0=linspace(-3,3,100)';
        x=zeros(length(x0),n);
        y=x;
        for j=1:n
            x(:,j)=x0+m0(1,j);
            cj=c0(j);
            y(:,j)=normpdf(x0,0,sqrt(cj));
        end
        plot(x,log10(y),':',x0+u,log10(normpdf(x0,0,k)*exp(g)),'k-');
        ylabel('Log10(pdf)');
    else
        if length(sc)<3
            if(sc(2)==1) & (n>1)    % covariance matrices are multiples of the identity
                sk=k*eye(d);
            else                    % diagonal covariance matrices
                sk=diag(k);
            end
            uk=eye(d);
            vk=uk;
        else                        % full covariance matrices
            [uk,sk,vk]=svd(k);
        end
        
        
        u2=uk(:,1:2);
        t0=linspace(0,2*pi,100);
        x=zeros(length(t0),n);
        y=x;
        x0=[cos(t0); sin(t0)];
        for j=1:n
            if length(sc)<3
                if(sc(2)==1) & (n>1)    % covariance matrices are multiples of the identity
                    cj=c0(j)*eye(2);
                else                    % diagonal covariance matrices
                    cj=u2'*diag(c0(:,j))*u2;
                end
            else                        % full covariance matrices
                cj=u2'*c0(:,:,j)*u2;
            end
            mj=u2'*m0(:,j);
            v=sqrt(sum((x0'/cj).*x0',2).^(-1));
            x(:,j)=mj(1)+v.*x0(1,:)';
            y(:,j)=mj(2)+v.*x0(2,:)';
        end
        
        if length(sc)<3
            if(sc(2)==1) & (n>1)    % covariance matrices are multiples of the identity
                cj=k*eye(2);
            else                    % diagonal covariance matrices
                cj=u2'*diag(k)*u2;
            end
        else                        % full covariance matrices
            cj=u2'*k*u2;
        end
        mj=u2'*u;
        v=sqrt(sum((x0'/cj).*x0',2).^(-1));
        plot(x,y,':',mj(1)+v.*x0(1,:)',mj(2)+v.*x0(2,:)','k-');
        axis equal;
    end
end
