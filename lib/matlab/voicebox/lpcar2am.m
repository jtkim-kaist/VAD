function [am,em]=lpcar2am(ar,p);
%LPCAR2AM Convert ar coefs to ar coef matrix [AM,EM]=(AR,P)
%AM is a 3-dimensional matrix of size (p+1,p+1,nf) where p is the lpc order
%and nf the number of frames.
%The matrix AM(:,:,*) is upper triangular with 1's on the main diagonal
%and contains the lpc coefficients for all orders from p down to 0.
%
%For lpc order p+1-r, AM(r,r:p+1,*), AM(p+1:-1:r,p+1,*) and EM(*,r) contain
%the lpc coefficients, reflection coefficients and the residual energy respectively.
%EM(:,1) is always 1.
%
%If A=am(:,:,*), R=toeplitz(rr(*,:)) and E=diag(em(*,:)), then
% A*R*A'=E; inv(R)=A'*(1/E)*A; A*R is lower triangular with the same diagonal as E
%
% For each j in 1:P we have AM(j:end:-1:j+1,*) = AM(j:end-1,end,*)'*am(j+1:end,j+1:end,*)
%
% Also em(*,1:end)' = em(*,2:end)'.*(1-am(1:end-1,end,*).^2)

% we should be able to do this more directly using the step down algorithm

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcar2am.m 713 2011-10-16 14:45:43Z dmb $
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

[nf,p1] = size(ar);
if any(ar(:,1)~=1)
  ar=ar./ar(:,ones(1,p1));
end
p0=p1-1;
if nargin < 2
   p=p0;
end
am=zeros(nf,p+1,p+1);
em=ones(nf,p+1);
e=ones(nf,1);
rf=ar;
if p>=p0
   for jj=1:p+1-p0
      am(:,jj:jj+p0,jj)=ar;
   end
else
   for j=p0:-1:p+2
      k = rf(:,j+1);
      d = (1-k.^2).^(-1);
      e = e.*d;
      wj=ones(1,j-1);
      rf(:,2:j) = (rf(:,2:j)-k(:,wj).*rf(:,j:-1:2)).*d(:,wj);
   end
   jj=0;
end
for jj=jj+1:p  
   j = p+2-jj;
   k = rf(:,j+1);
   d = (1-k.^2).^(-1);
   e = e.*d;
   wj=ones(1,j-1);
   rf(:,2:j) = (rf(:,2:j)-k(:,wj).*rf(:,j:-1:2)).*d(:,wj);
   am(:,jj:end,jj)=rf(:,1:j);
   em(:,jj)=e;
end
em(:,end)=e./(1-rf(:,2).^2);
am(:,end,end)=1;
am=permute(am,[3 2 1]);
   
