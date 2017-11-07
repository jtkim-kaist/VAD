function [ar,arp,aru,g]=lpcrf2ar(rf)
%LPCRF2AR Convert reflection coefs to autoregressive coefs [AR,ARP,ARU,G]=(RF)
%
% Input:  RF(:,p+1) gives reflection coefficients of one or more p-section lossless tubes 
% Ouputs: G is the gain of the all-pole AR filter
%         AR/G is the transfer function from U_in to the glottal input wave, U_g.
%               AR(:,1)=1 always.
%         ARP*K is the transfer function from U_in to the pressure just after the glottis
%               where K = rho*c/Alips: rho = air density 1.23 kg/m^3, c=sound speed 340 m/s, 
%               Alips = effective area of free space beyond the lips.
%         ARU is the transfer function from U_in to the total volume velocity through the glottis
% 
%              where U_in=z^(p/2)*U_lips is the time-advanced volume velocity at the lips
%
%         Energy into the vcal tract is equal to K*filter(ARP,1,Ulips).*filter(ARU,1,Ulips)
%              reverse glottal flows divided by 1-r0 where r0 is the glottal reflection coefficient.
%              The scale factor is included to avoid a zero answer when the glottis is closed giving r0=1.
%
% The transfer functions have ar(:,1)=art(:,1)=1
% They should both be multiplied by z^(p/2)/prod(1+rf) to correct the absolute
% gain and to compensate for the delay of p/2 samples along the length of the tube.
%
% The energy into the vocal tract is given by ars(speech) * are(speech)
%
% Ref: D. M. Brookes and H. P. Loke. "Modelling energy flow in the vocal tract with
%           applications to glottal closure and opening detection." In Proc ICASSP'99, pages 213-216, Mar 1999.


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcrf2ar.m 713 2011-10-16 14:45:43Z dmb $
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
p2=p1+1;
p=p1-1;
pm=p-1;
arf=[ones(nf,1) zeros(nf,p)];
arr=[zeros(nf,p) rf(:,p1)];
cr=zeros(nf,p);
for k=1:p-1
  rk=rf(:,(p1-k)*ones(1,k));
  cr(:,1:k)=arr(:,p2-k:p1);
  arr(:,p1-k:p)=arr(:,p1-k:p)+rk.*arf(:,1:k);
  arf(:,2:k+1)=arf(:,2:k+1)+rk.*cr(:,1:k);
end
r1=rf(:,1);
ar=arf+r1(:,ones(1,p1)).*arr;
if nargout>1
   kp=prod(1-rf(:,2:p1),2);
   arp=(arf-arr)./kp(:,ones(1,p1));
   if nargout>2
      g=prod(1+rf(:,2:p1),2);
      aru=(arf+arr)./g(:,ones(1,p1));
      if nargout>3
         g=g.*(1+rf(:,1));
      end
   end
end


