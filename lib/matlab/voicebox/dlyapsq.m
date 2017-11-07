function v=dlyapsq(a,b)
%DLYAPSQ Solves the discrete Lyapunov equation AV'VA' - V'V + BB' = 0
% V is upper triangular with real non-negative diagonal entries
% this is equivalent to v=chol(dlyap(a,b*b')) but better conditioned numerically

%      Copyright (C) Mike Brookes 2002
%      Version: $Id: dlyapsq.m 7339 2016-01-06 18:05:30Z dmb $
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

[q,s]=schur(a');
[q,s]=rsf2csf(q,s);
[qd,r]=qr(b'*q,0);
% save r for testing
r0=r;
[m,n]=size(r);
u=zeros(n,n);
if m==1
   for i=1:n-1
      in=i+1:n;
      si=s(i,i);
      aa=sqrt(1-si*si');
      u(i,i)=r(1)/aa;
      u(i,in)=(u(i,i)*si'*s(i,in)+aa*r(2:end))/(eye(n-i)-si'*s(in,in));
      r=aa*(u(i,i)*s(i,in)+u(i,in)*s(in,in))-si*r(2:end);
   end
   u(n,n)=r/sqrt(1-s(n,n)*s(n,n)');
   
else
   w=zeros(m,1); w(m)=1;
   em=eye(m);
   for i=1:n-m
      in=i+1:n;
      si=s(i,i);
      aa=sqrt(1-si*si');
      u(i,i)=r(1,1)/aa;
      u(i,in)=(u(i,i)*si'*s(i,in)+aa*r(1,2:end))/(eye(n-i)-si'*s(in,in));
      vv=aa*(u(i,i)*s(i,in)+u(i,in)*s(in,in))-si*r(1,2:end);
      rr=zeros(m,n-i);
      rr(1:m-1,:)=r(2:end,2:end);
      [qq,r]=qrupdate(em,rr,w,vv');
   end
   for i=n-m+1:n-1
      in=i+1:n;
      si=s(i,i);
      aa=sqrt(1-si*si');
      u(i,i)=r(1,1)/aa;
      u(i,in)=(u(i,i)*si'*s(i,in)+aa*r(1,2:end))/(eye(n-i)-si'*s(in,in));
      vv=aa*(u(i,i)*s(i,in)+u(i,in)*s(in,in))-si*r(1,2:end);
      rr=zeros(n-i+1,n-i);
      rr(1:n-i,:)=r(2:end,2:end);
      [qq,rr]=qrupdate(eye(n-i+1),rr,w(m-n+i:end),vv');
      r=rr(1:n-i,:);
   end
   
   u(n,n)=r/sqrt(1-s(n,n)*s(n,n)');
   
end

v=triu(qr(u*q'));
dv=diag(v);
ix=dv~=0;
v(ix,:)=diag(abs(dv(ix))./dv(ix))*v(ix,:);
if isreal(a) & isreal(b)
   v=real(v);
end
