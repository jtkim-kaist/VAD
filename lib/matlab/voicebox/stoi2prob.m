function p=stoi2prob(s,m)
%STOI2PROB convert STOI to probability
%
%  Inputs:  S(M,N)       matrix containing STOI values
%           M            mapping: 'i' IEEE sentences [default]
%                                 'd' Dantale corpus
%
% Outputs:  P(M,N)       Corresponding probability values
%
% STOI is an intelligibility metric described in [1]. The
% mapping from STOI to intelligibilty is corpus-dependent:
% this functions implements two mappings given in [1].
%
% [1]	C. H. Taal, R. C. Hendriks, R. Heusdens, and J. Jensen.
%       An algorithm for intelligibility prediction of time–frequency weighted noisy speech.
%       IEEE Trans. Audio, Speech, Language Processing, 19 (7): 2125–2136, 2011.
%       doi: 10.1109/TASL.2011.2114881.

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: stoi2prob.m 5887 2015-03-18 09:42:18Z dmb $
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

if nargin<2 || ~numel(m) || m(1)=='i'
    a=-17.4906;
    b=9.6921;
    if ~nargin
        s=[];
    end
else
    a=-14.5435;
    b=7.0792;
end
p=1./(1+exp(a*s+b));
if ~nargout
    if ~nargin
        s=linspace(0,1,100);
    end
    plot(s,stoi2prob(s,'d'),':k',s,stoi2prob(s,'i'),'-k')
    xlabel('STOI');
    ylabel('Intelligibility');
    legend('Dantale corpus','IEEE sentences','location','best');
end