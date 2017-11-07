function [cc,c0]=lpcpf2cc(pf,np,f)
%LPCPF2CC Convert power spectrum to complex cepstrum CC=(PF,NP)
%
%  Inputs: pf(nf,n)    Power spectrum, uniformly spaced DC to Nyquist
%          np          Number of cepstral coefficients to calculate [n-1]
%          f(1,n)      Frequencies of pf columns [linspace(0,0.5,n)]
%
% Outputs: cc(nf,np)  Complex spectrum from DC to Nyquist
%          c0(nf,1)   The zeroth cepstral coefficient, c(0)
%
% Note since the log spectrum is not normally bandlimited, this conversion is not exact unless n >> np

%      Copyright (C) Mike Brookes 1998-2014
%      Version: $Id: lpcpf2cc.m 5025 2014-08-22 17:07:24Z dmb $
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
[nf,nq]=size(pf);
if nargin<2 || ~numel(np) np=nq-1; end
if nargin<3 || ~numel(f)
    cc=rsfft(log(pf.')).'/(2*nq-2);
    c0=cc(:,1)*0.5;
    cc(:,nq)=cc(:,nq)*0.5;
    if np>nq-1
        cc=[cc(:,2:nq) zeros(nf,np-nq+1)];
    else
        cc=cc(:,2:np+1);
    end
else
    cc=0.5*(log(pf)/cos(2*pi*f(:)*(0:min(np,nq-1))));
    c0=cc(:,1);
    cc=cc(:,2:end);
    if np>nq-1
        cc(1,np)=0;
    end
end


