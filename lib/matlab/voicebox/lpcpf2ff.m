function [ff,fo]=lpcpf2ff(pf,np,fi)
%LPCPF2FF Convert power spectrum to complex spectrum [FF,FO]=(PF,NP,FI)
%
%  Inputs: pf(nf,n)     Power spectrum at n discrete frequencies, one frame per row
%          np           Number of complex cepstral coefficients to use (excluding c0) [n-1]
%                          should be greater than the sum of the numerator
%                          and denominator filter orders but less than n
%          fi(1,n)      Vector of frequencies [linspace(0,0.5,n)]
%                         including this argument slows down the routine
%
% Outputs: ff(nf,n)     Complex spectrum (pf = abs(ff).^2
%          fo(1,n)      Vector of frequencies
%
% This routine converts a power spectrum into the corresponding complex
% spectrum. It determines the phase spectrum under the assumption that it
% is minimum phase. The routine works by converting first to the compex
% cepstrum.

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: lpcpf2ff.m 5026 2014-08-22 17:47:43Z dmb $
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
if nargin<3 fi=[];
    if nargin<2
        np=nq-1; % number of cepstal coefficients (excl c(0))
    end
end
[cc,c0]=lpcpf2cc(pf,np,fi);
if ~numel(fi)
    fi=nq-1;
end
[fx,fo]=lpccc2ff(cc,fi,-1,c0);
ff=sqrt(pf).*exp(1i*angle(fx));
if ~nargout
    subplot(2,1,2);
    plot(fo,unwrap(angle(ff.')));
    xlabel('Normalized frequency f/f_s');
    ylabel('Phase (rad)');
    subplot(2,1,1);
    plot(fo,db(abs(ff.')),'-b',fo,db(pf.')/2,':k');
    xlabel('Normalized frequency f/f_s');
    ylabel('Gain (dB)');
end


