function [a,d,e]=dualdiag(w,b)
%DUALDIAG Simultaneous diagonalisation of two hermitian matrices [A,D,E]=(W,B)
% Inputs:   W,B     Two square hermitian matrices
%
% Outputs:  A       Diagonalizing matrix (not normally unitary or hermitian)
%           D       Real diagonal matrix elements: A'*B*A=diag(D) (see note below)
%           E       Real diagonal matrix elements: A'*W*A=diag(E)
%
% Note: At least one of W and B must be either positive definite or negative
% definite. If this is not the case, then D=A'*B*A and E=A'*W*A will be
% complex hermitian matrices rather than being vectors of real diagonal elements.
%
% The columns of A will be ordered so that abs(D./E) is in descending order.
% If two output arguments are given then A will be scaled to make diag(E)=I
% but, if W is singular, this will cause some elements of A to be infinite.
%
% If A'*B*A=diag(D) and A'*W*A=diag(E) then A'*W*A*diag(1./E)=I so A'*B*A=A'*W*A*diag(D./E)
% and hence B*A=W*A*diag(D./E) so the columns of A are the eigenvectors of W\A or
% equivalently the generalized eigenvectors of (B,W).
%
% Suppose we have several N-dimensional data row-vectors arising from each of C different classes of data.
% for each class, c, we can form the mean data vector m(c) and the within-class covariance matrix W(c)
% We can then form the between class covariance matrix B by taking the covariance of the mean vectors m(1), m(2), ...
% and also the averaged within-class covariance matrix W by averaging W(1), W(2), ...
% If we then take A=dualdiag(W,B) and postmultiply all our original data vectors by A, we obtain new
% data vectors for which the average within-class covariance matrix is the identity and for which
% the first few components contain most of the information that is useful in discriminating between classes.

%      Copyright (C) Mike Brookes 1997-2013
%      Version: $Id: dualdiag.m 2867 2013-04-02 11:27:18Z dmb $
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
[a,l]=eig(b+b',w+w');               % generalized eigendecomposition
if isreal(l)
    [d,i]=sort(abs(diag(l)),'descend'); % sort by absolute value
    if nargout==2
        q=sqrt(diag(a'*w*a))'.^(-1);    % scale to make e=1
        a=a(:,i).*q(ones(size(w,1),1),i); % reorder and scale columns of a
    else
        a=a(:,i);                       % reorder columns of a to match eigenvalues
        e=real(diag(a'*w*a));
    end
    d=real(diag(a'*b*a));
else
    d=a'*b*a;
    e=a'*w*a;
end
