function q=importsii(f,m)
%IMPORTSII calculates the SII importance function per Hz or per Bark Q=(F,M)
% Inputs:
%        f(n)   Frequencies to which to calculate importance in Hz
%               or Bark according to 'b' flag.
%        m      Mode string with some of the following flags:
%                'b'  Frequencies given in Bark rather than Hz
%                'c'  Calculate cumulative importance for f<f(i)
%                'd'  Calculate importance of n-1 bands: band i is f(i) to f(i+1)
%                'h'  Calculate importance per Hz or per Bark (accoring to 'b' flag)
% Outputs:
%        q(n) or q(n-1) gives the importance at each of the f(i) or else,
%               if the 'd' flag is specified, in the band from f(i) to f(i+1).

% The importance function is based on the piecewise linear function
% defined in Fig 3 of [2]. This is integrated to give the cumulative
% importance function. It is modified slightly from Fig 3 so that the
% constant portion extends from 4 to 18 Bark (critical bands 5 to 18).
% we then fit a linear portion at either end to force the correct integral
% and ensure continuity of the importance at 4 and 18. The importance
% function is zero outside the range [1.286 21.948] bark or [130.1 9361]
% Hz.
%
% References:
%  [1]  Methods for the calculation of the speech intelligibility index.
%       ANSI Standard S3.5-1997 (R2007), American National Standards Institute, 1997.
%  [2]  C. V. Pavlovic. Derivation of primary parameters and procedures for use in
%       speech intelligibility predictions. J. Acoust Soc Amer, 82: 413–422, 1987.

%	   Copyright (C) Mike Brookes 2006
%      Version: $Id: importsii.m 5099 2014-09-09 14:39:17Z dmb $
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

persistent mi ci ai bi xi0 xi1
if isempty(mi)
    % cumulative importance function =  mi*b+c+ai*(b<4)*(b-4)^2-bi*(b>18)*(b-18)^2
    % for x in (xi0,xi1)
    ci4=0.0783;
    ci18=0.8861;
    mi=(ci18-ci4)/14;
    ci=ci4-4*mi;
    ai=mi^2/(4*(4*mi+ci));
    bi=mi^2/(4*(1-18*mi-ci));
    xi0=4-mi/(2*ai);
    xi1=18+mi/(2*bi);
end
if nargin<2
    m=' ';
end
if any(m=='b')
    b=f;
else
    [b,d]=frq2bark(f);
end
if any(m=='c') || any(m=='d')
    q=mi*b+ci+ai*(b<4).*(b-4).^2-bi*(b>18).*(b-18).^2;
    q(b<xi0)=0;
    q(b>xi1)=1;
    if any(m=='d')
        q=q(2:end)-q(1:end-1);
    end
else
    q=mi+ai*(b<4).*(b-4)-bi*(b>18).*(b-18);
    q(b<xi0)=0;
    q(b>xi1)=0;
    if ~any(m=='b')
        q=q./d;
    end
end
if ~nargout
    if any(m=='d')
        ix=(1:2*length(q))/2;
        plot(f(1+floor(ix)),q(ceil(ix)));
    else
        plot(f,q);
    end
    if any(m=='b')
        xlabel('Frequency (Bark)');
    else
        xlabel('Frequency (Hz)');
    end
    ylabel('Importance');
    if any(m=='c')
        title('SII Cumulative Importance');
    elseif any(m=='d')
        title('SII Band Importance');
    else
        title('SII Importance Function');
    end
end

