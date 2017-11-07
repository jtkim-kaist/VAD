function [a,f,B,W]=ldatrace(b,w,n,c)
%LDATRACE Calculates an LDA transform to maximize trace discriminant [a,f,B,W]=(b,w,n,c)
% If a feature vector X can come from one of several class and W and B are respectively
% the within-class and between-class covariance matrices, then the generalized Fisher discriminant
% F=trace(W\B) is a measure of how well the feature vector discriminates between the classes.
% If we choose a rectangular (tall, skinny) transformation matrix, we can define a smaller
% feature vector Y=A'*X. The aim of this routine is to choose A to maximize the Fisher
% discriminant. We assume that W is positive definite and B is positive semi-definite.
% The input argument C allows the uset to pre-specify some of the columns of A.
%
% Inputs:
%     w[m,m] = within class covariance matrix of x
%     b[m,m] = between class covariance matrix of x [default = I]
%     n is the number of columns in output matrix A [default = M]
%     c[m,r] specifies the first few columns of A to be predefined values [default = null)
%
% Outputs:
%     a[m,n] is the transformation matrix: y=a'*x
%     f[n,1] gives the incremental gain in f value for successive columns of A 
%     B(n,n) gives the between-class covariance matrix of y
%     W[n,n] gives the within-class covariance matrix of y 

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: ldatrace.m 713 2011-10-16 14:45:43Z dmb $
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

m=size(b,1);    % dimension of data vectors
if nargin<4
    r=0;
    if nargin<3
        n=m;
        if nargin<2
            w=eye(m);
        end
    end
else
    r=size(c,2);    % number of columns that are pre-specified
end
if r
    if n>r          % need to find additional vectors
        g=chol(w);
        v=g\null(c'*g');
        [p,l,q]=svd(v'*b*v);
        a(:,r+1:n)=v*p(:,1:n-r);
        a(:,1:r)=c;
    else
        a=c;        % no new vectors to find
    end
    if nargout>1
        ari=a/triu(qr(chol(a'*w*a))); % matrix a must be of full rank
        f=diag(ari'*b*ari);
    end
else
    [g,d]=eig(b,w,'qz');
    [ds,is]=sort(-diag(d));
    a=g(:,is(1:n));
    if nargout>1
        f=-ds(1:n);
    end
end
if nargout > 2
    B=a'*b*a;
    if nargout > 3
        W=a'*w*a;
    end
    
end
