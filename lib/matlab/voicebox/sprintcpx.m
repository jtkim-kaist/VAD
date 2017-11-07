function s=sprintcpx(z,f)
%SPRINTCPX  format a complex number for printing S=(Z,F)
%
% Usage: fprintf('%s',sprintcpx(z));
%
%  Inputs: z   a complex number to print
%          f   optional formatting string as in fprintf e.g. '0.2f' [default: 'g']
%              may also include 'i' or 'j' [default] to control sqrt(-1) symbol.
%
% Outputs: s   formatted output string

%      Copyright (C) Mike Brookes 2015
%      Version: $Id: sprintcpx.m 7297 2015-12-08 12:27:09Z dmb $
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

if nargin<2 || ~numel(f)
    f='g';
end
if any(f=='i')
    ij='i';
else
    ij='j';
end
f((f=='i')|(f=='j'))=[]; % remove i and j specifiers
if ~numel(f)
    f='g';
end
if any(f=='+')
    pl='';
else
    pl='+';
end
f=['%' f];
a=real(z);
b=imag(z);
jx=[1 3 2 4 3 4 1 3 2];
ix=jx(3*sign(a)+sign(b)+5);
switch(ix)
    case 1
        s=sprintf([f f ij],a,b);
    case 2
        s=sprintf([f pl f ij],a,b);
    case 3
        s=sprintf(f,a);
    case 4
        s=sprintf([f ij],b);
end
