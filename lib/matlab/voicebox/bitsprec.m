function y=bitsprec(x,n,mode)
%BITSPREC round values to a specified fixed or floating precision (X,N,MODE)
%
% mode is of the form 'uvw' where:
%     u: s - n significant bits (default) 
%        f - fixed point: n bits after binary point
%     v: n - round to nearest (default)
%        p - round towards +infinity
%        m - round towards -infinity
%        z - round towards zero
% w is only needed if v=n in which case it dictates what to
%        do if x is min-way between two rounded values:
%     w: p,m - as above
%        e - round to nearest even number (default)
%        o - round to nearest odd number
%        a - round away from zero
% mode='*ne' and '*no' are convergent rounding and introduce
% no DC offset into the result so long as even and odd integer parts are
% equally common.
%
% Examples of y=bitsprec(x,0,'***'):
%
%      x    fp-  fm-  fz-  fne  fno  fnp  fnm  fna 
%   
%     2.5    3    2    2    2    3    3    2    3
%     1.5    2    1    1    2    1    2    1    2
%     1.1    2    1    1    1    1    1    1    1
%     1.0    1    1    1    1    1    1    1    1
%     0.9    1    0    0    1    1    1    1    1
%     0.5    1    0    0    0    1    1    0    1
%     0.1    1    0    0    0    0    0    0    0
%    -0.1    0   -1    0    0    0    0    0    0
%    -0.5    0   -1    0    0   -1    0   -1   -1
%    -0.9    0   -1    0   -1   -1   -1   -1   -1
%    -1.5   -1   -2   -1   -2   -1   -1   -2   -2

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: bitsprec.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<3
   mode='sne';
end
if mode(1)=='f'
   e=0;
else
   [x,e]=log2(x);
end
switch mode(2)
case 'p'
   y=pow2(ceil(pow2(x,n)),e-n);
case 'm'
   y=pow2(floor(pow2(x,n)),e-n);
case 'z'
   y=pow2(fix(pow2(x,n)),e-n);
otherwise
   switch mode(3)
   case 'a'
      y=pow2(round(pow2(x,n)),e-n);
   case 'p'
      y=pow2(floor(pow2(x,n)+0.5),e-n);
   case 'm'
      y=pow2(ceil(pow2(x,n)-0.5),e-n);
   otherwise
      z=pow2(x,n-1);
      switch mode(3)
      case 'e'
         y=pow2(floor(pow2(x,n)+0.5)-floor(z+0.75)+ceil(z-0.25),e-n);
      case 'o'
         y=pow2(ceil(pow2(x,n)-0.5)+floor(z+0.75)-ceil(z-0.25),e-n);
      end      
   end
end

