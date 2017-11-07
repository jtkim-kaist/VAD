function s=phon2sone(p)
%PHON2SONE convert PHON loudness values to SONEs s=(p)
%Inputs:    p is a matrix of phon values
%
%Outputs:   s is a matrix, the same size as p, of sone values
%
% The phon scale measures perceived loudness in dB; at 1 kHz it is identical to dB SPL
% relative to 20e-6 Pa sound pressure. The sone scale is proportional to apparent loudness
% and, by definition, equals 1 at 40 phon. The form of the loudness curve is taken from [1].
% The hearing threshold at 1 kHz for 18 to 25 year olds with normal hearing is taken from [2].
%
% Refs: [1]	J. Lochner and J. Burger. Form of the loudness function in the presence of masking noise.
%           The Journal of the Acoustical Society of America, 33: 1705, 1961.
%       [2]	ISO/TC43. Acoustics – normal equal-loudness-level contours.
%           Standard ISO 226:2003, Aug. 2003.


%      Copyright (C) Mike Brookes 2012-2013
%      Version: $Id: phon2sone.m 3295 2013-08-02 14:03:11Z dmb $
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
persistent a b d
if isempty(a)
    b=log(10)*0.1*0.27; % 0.27 is the exponent from [1] and [2]
    d=exp(b*2.4); % 2.4 dB is teh hearing threshold from [2]
    a=1/(exp(b*40)-d); % scale factor to make p=40 give s=1
end
if nargout>0
    s=a*(exp(b*p)-d);
else
    if nargin<1 || isempty(p)
        pp=linspace(5,90,100)'; % phon values
    else
        pp=p;
    end
    semilogy(pp,phon2sone(pp));
    axisenlarge(-1);
    yticksi;
    xlabel('phon = dB SPL @ 1 kHz');
    ylabel('sone');
end
