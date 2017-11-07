function [n,f,a,b]=lpcar2fm(ar,t)
%LPCAR2RF Convert autoregressive coefficients to formant freq+amp+bw [N,F,A,B]=(AR,T)
%
% Input:   ar(:,p+1)  Autoregressive coefficients
%          t          Threshold (see below)
% Output:  n          Number of formants found
%          f          Formant frequencies in normalized Hz (in increasing order)
%          a          Formant amplitudes
%          b          Formant bandwidths in normalized Hz
%
% The number of columns in the output arrays f, a and b is max(n); surplus positions
% in any given row have f=b=0.
%
% In determining formants, poles are ignored if any of the following hold:
%        (a) they are on the real axis
%        (b) they have bandwidth > t*frequency (if t>0)
%        (c) they have bandwidth > -t (if t<=0)

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcar2fm.m 713 2011-10-16 14:45:43Z dmb $
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

[nf,p1]=size(ar);
p=p1-1;
d=(1:nf)';
zz=lpcar2zz(ar);
ig=imag(zz)<=0;
n=p1-1-sum(ig,2);
mn=max(n);

% remove redundant columns

if mn<p
   [ig,ix]=sort(ig,2);
   zz=reshape(zz(d(:,ones(1,mn))+nf*(ix(:,1:mn)-1)),nf,mn);
   ig(:,mn+1:end)=[];
end

zz(ig)=1;      % to prevent infinities
f=angle(zz)*0.5/pi;
b=-log(abs(zz))/pi;
if nargin > 1
   if t>0
      ig=ig | b>t*f;
   else
      ig=ig | b+t>0;
   end
end
f(ig)=0;
b(ig)=0;
n=mn-sum(ig,2);
m=max(n);

% remove redundant columns

[igf,ix]=sort(ig+f,2);
dd=d(:,ones(1,m))+nf*(ix(:,1:m)-1);
zz=reshape(zz(dd),nf,m);
f=reshape(f(dd),nf,m);
b=reshape(b(dd),nf,m);
ig=reshape(ig(dd),nf,m);

% now calculate gain
ap=permute(ar,[1 3 2]);
pw=permute(-2*pi*1i*(0:p),[1 3 2]);
a=abs(sum(ap(:,ones(1,m),:).*exp(pw(ones(1,nf),ones(1,m),:).*f(:,:,ones(1,p1))),3)).^(-1);
