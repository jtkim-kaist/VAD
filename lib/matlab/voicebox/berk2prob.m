function [p,d]=berk2prob(b)
%BERK2PROB convert Berksons to probability
%
%  Inputs:  B(M,N)       matrix containing Berkson values
%
% Outputs:  P(M,N)       Corresponding probability values
%           D(M,N)       Corresponding derivatives dP/dB
%
% Berksons, or log-odds, are a nonlinear scale for measuring
% probability defined by B = log2(P./(1-P)).
% When Berksons are used to measure probability, a logistic
% psychometric function becomes linear.
%
% The inverse function is berk2prob()

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: berk2prob.m 4501 2014-04-24 06:28:21Z dmb $
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
p=1-1./(1+pow2(b));
if nargout>1
    d=log(2)*p.*(1-p);
end