function [rr,ar]=lpcrf2rr(rf,p);
%LPCRR2AR convert reflection coefs to autocorrelation coefs [RR,AR]=(RF,P)
%
% Inputs:  rf(:,n+1)  reflection coefficients: one row per frame
%          p          specifies number of rr coefficients to calculate (default=n)
% Outputs: rr(:,p+1)  autocorrelation coefficients
%          ar(:,n+1)  AR filter coefficients

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcrf2rr.m 2460 2012-10-29 22:20:45Z dmb $
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

[nf,p1]=size(rf);
p0=p1-1;
if p0
   a = rf(:,2);
   rr=[ones(nf,1) -a zeros(nf,p0-1)];
   e = (a.^2-1);
   for n = 2:p0
      k=rf(:,n+1);
      rr(:,n+1) =k.*e - sum(rr(:,n:-1:2).*a,2);
      a = [a+k(:,ones(1,n-1)).*a(:,n-1:-1:1) k];
      e = e.*(1-k.^2);
   end
   ar = [ones(nf,1) a];
   r0=sum(rr.*ar,2).^(-1);
   rr=rr.*r0(:,ones(1,p1));
   if nargin>1 && ~isempty(p)
      if p<p0
         rr(:,p+2:p1)=[];
      else
         rr=[rr zeros(nf,p-p0)];
         af=-ar(:,p1:-1:2);
         for i=p0+1:p
            rr(:,i+1)=sum(af.*rr(:,i-p0+1:i),2);
         end
      end
   end
else
   rr=ones(nf,1);
   ar=rr;
end

