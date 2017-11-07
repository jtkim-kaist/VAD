function [am,em]=lpcrr2am(rr);
%LPCRR2AM Convert autocorrelation coefs to ar coef matrix [AM,EM]=(RR)
%AM is a 3-dimensional matrix of size (p+1,p+1,nf) where p is the lpc order
%and nf the number of frames.
%The matrix AM(:,:,*) is upper triangular with 1's on the main diagonal
%and contains the lpc coefficients for all orders from p down to 0.
%
%For lpc order p+1-r, AM(r,r:p+1,*), AM(p+1:-1:r,p+1,*) and EM(*,r) contain
%the lpc coefficients, reflection coefficients and the residual energy respectively.
%
%If A=am(:,:,*), R=toeplitz(rr(*,:)) and E=diag(em(*,:)), then
% A*R*A'=E; inv(R)=A'*(1/E)*A; A*R is lower triangular with the same diagonal as E
%
% This routine is equivalent to: c=chol(inv(toeplitz(rr))); d=diag(c).^-1; em=d.^2; am=diag(d)*c

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcrr2am.m 713 2011-10-16 14:45:43Z dmb $
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

[nf,p1]=size(rr);
p=p1-1;
p2=p1+1;
am=zeros(nf,p1,p1);
em=zeros(nf,p1);
am(:,p1,p1)=1;
em(:,p1)=rr(:,1);
ar=ones(nf,p1);
ar(:,2) = -rr(:,2)./rr(:,1);
e = rr(:,1).*(ar(:,2).^2-1);
for n = 2:p
   q=p2-n;
   em(:,q)=-e;
   am(:,q:p1,q)=ar(:,1:n);
   k = (rr(:,n+1)+sum(rr(:,n:-1:2).*ar(:,2:n),2)) ./ e;
   ar(:,2:n) = ar(:,2:n)+k(:,ones(1,n-1)).*ar(:,n:-1:2);
   ar(:,n+1) = k;
   e = e.*(1-k.^2);
end
em(:,1)=-e;
am(:,:,1)=ar;
am=permute(am,[3 2 1]);

