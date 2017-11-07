function p=mintrace(x)
%MINTRACE find row permutation to minimize the trace p=(x)
%Inputs:    x(n,n)  is a square real-valued matrix
%
%Outputs:   p(n,1)  is the row permutation that minimizes the trace
%                   so that trace(x(p,:)) is as small as possible

%      Copyright (C) Mike Brookes 2012-2013
%      Version: $Id: mintrace.m 2718 2013-02-23 09:34:45Z dmb $
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
n=size(x,1);    % input x must be square
p=permutes(n);  % try all permutations
c=0:n:n^2-1;    % convert olumns to single indexing
[v,i]=min(sum(x(p+c(ones(length(p),1),:)),2));
p=p(i,:)';

