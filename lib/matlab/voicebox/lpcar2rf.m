function rf=lpcar2rf(ar)
%LPCAR2RF Convert autoregressive coefficients to reflection coefficients AR=(RF)
%
% Input:   ar(:,p+1)  Autoregressive coefficients
% Output:  rf(:,p+1)  Reflection coefficients with rf(:,1)=1


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcar2rf.m 8211 2016-07-20 20:59:16Z dmb $
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
if p1==1
   rf=ones(nf,1);
else
   if any(ar(:,1)~=1)
      ar=ar./ar(:,ones(1,p1));
   end
   rf = ar;
   for j = p1-1:-1:2
      k = rf(:,j+1);
      d = (1-k.^2).^(-1);
      wj=ones(1,j-1);
      rf(:,2:j) = (rf(:,2:j)-k(:,wj).*rf(:,j:-1:2)).*d(:,wj);
   end
end

