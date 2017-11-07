function [m,a]=lpcstable(ar)
%LPCSTABLE Test AR coefficients for stability and stabilize if necessary [MA,A]=(AR)
%
% Usage: (1) [m,ar]=lpcstable(ar); % force ar polynolials to be stable
%
%   Input:  ar(:,p+1)  Autoregressive coefficients
% Outputs:  m(:,1)    Mask identifying stable polynomials
%           a(:,p+1)  Stabilized polynomials formed by reflecting unstable
%                       poles in unit circle (with a(:,1)=1)

%      Copyright (C) Mike Brookes 2016
%      Version: $Id: lpcstable.m 8558 2016-09-22 08:22:54Z dmb $
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
mm=ar(:,1)~=1;
if any(mm)          % first ensure leading coefficient is always 1
    ar(mm,:)=ar(mm,:)./ar(mm,ones(1,p1));
end
if p1==1
    m=ones(nf,1); % 0'th order filter is always stable
elseif p1==2
    m=abs(ar(:,2))<1; % 1'st order filter
else
    rf = ar;
    k = rf(:,p1); % check final coefficient in range
    m=abs(k)<1;
    if any(m)
        d = (1-k(m).^2).^(-1);
        wj=ones(1,p1-2);
        rf(m,2:p1-1) = (rf(m,2:p1-1)-k(m,wj).*rf(m,p1-1:-1:2)).*d(:,wj);
        for j = p1-2:-1:2
            k(m) = rf(m,j+1);
            m(m)=abs(k)<1;
            if ~any(m), break, end
            d = (1-k(m).^2).^(-1);
            wj=ones(1,j-1);
            rf(m,2:j) = (rf(m,2:j)-k(m,wj).*rf(m,j:-1:2)).*d(:,wj);
        end
    end
end
if nargout>1
    a=ar;
    if ~all(m)
        for i=find(~m)'                 % unstable frames
            z=roots(a(i,:));
            k=abs(z)>1;                 % find any unstable roots
            z(k)=conj(z(k)).^(-1);      % invert them
            a(i,:)=real(poly(z));       % force a real polynomial
        end
    end
end
