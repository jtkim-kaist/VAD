function [ar,e,k]=lpcauto(s,p,t)
%LPCAUTO  performs autocorrelation LPC analysis [AR,E,K]=(S,P,T)
%  Inputs:
%     s(ns)   is the input signal
%	   p       is the order (default: 12)
%	   t(nf,3) specifies the frames size details: each row specifies
%	           up to three values per frame: [len anal skip] where:
%		           len     is the length of the frame (default: length(s))
%		           anal    is the analysis length (default: len)
%		           skip    is the number of samples to skip at the beginning (default: 0)
%	           If t contains only one row, it will be used repeatedly
%	           until there are no more samples left in s.
%
% Outputs:
%          ar(nf,p+1) are the AR coefficients with ar(1) = 1
%          e(nf)      is the energy in the residual.
%                     sqrt(e) is often called the 'gain' of the filter.
%          k(nf,2)    gives the first and last sample of the analysis interval

% Notes:
%
% (1) The first frame always starts at sample s(1) and the analysis window starts at s(t(1,3)+1).
% (2) The elements of t need not be integers.
% (3) The analysis interval is always multiplied by a hamming window
% (4) As an example, if p=3 and t=[10 5 2], then the illustration below shows
%     successive frames labelled a, b, c, ... with capitals for the
%     analysis regions. Note that the first frame starts at s(1)
%
%	  a a A A A A A a a a b b B B B B B b b b c c C C C C C c c c d ...
%
%     For speech processing, it can be advantageous to restrict the analysis regions
%     to time intervals when the glottis is closed.
%
% (5) Frames can overlap: e.g. t=[10 20] will use analysis frames of
%     length 20 overlapped by 10 samples.
% (6) For speech processing p should be at least 2*f*l/c where f is the sampling
%     frequency, l the vocal tract length and c the speed of sound. For a typical
%     male (l=17 cm) this gives f/1000.

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: lpcauto.m 713 2011-10-16 14:45:43Z dmb $
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

s = s(:); % make it a column vector
if nargin < 2 p=12; end;
if nargin < 3 t=length(s); end;
%if nargin < 4 w='ham'; end;
[nf,ng]=size(t);
if ng<2 t=[t t]; end;
if ng<3 t=[t zeros(nf,1)]; end;
if nf==1
    nf=floor(1+(length(s)-t(2)-t(3))/t(1));
    tr=0;
else
    tr=1;
end;

ar=zeros(nf,p+1);
ar(:,1)=1;
e=zeros(nf,1);

t1=1;
it=1;
nw=-1;
zp=zeros(1,p);
r=(0:p);
for jf=1:nf
    k(jf,1) = ceil(t1+t(it,3));
    k(jf,2) = ceil(t1+t(it,3)+t(it,2)-1);
    cs = (k(jf,1):k(jf,2)).';
    nc = length(cs);
    pp=min(p,nc);
    dd=s(cs);
    if nc~=nw
        % possibly we should have a window whose square integral equals unity
        ww=hamming(nc); nw=nc;
        y=zeros(1,nc+p);
        c=(1:nc)';
    end
    wd=dd(:).*ww;        % windowed data vector
    y(1:nc)=wd;          % data vector with p appended zeros
    z=zeros(nc,pp+1);    % data matrix
    %  was previously  z(:)=y(c(:,ones(1,pp+1))+r(ones(nc,1),1:pp+1));
    z(:)=y(repmat(c,1,pp+1)+repmat(r,nc,1));
    rr=wd'*z;
    rm=toeplitz(rr(1:pp));
    rk=rank(rm);
    if rk
        if rk<pp
            rm=rm(1:rk,1:rk);
        end
        ar(jf,2:rk+1)=-rr(2:rk+1)/rm;
    end
    e(jf)=rr*ar(jf,1:pp+1)';
    t1=t1+t(it,1);
    it=it+tr;
end

