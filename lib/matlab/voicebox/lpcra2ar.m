function ar=lpcra2ar(ra,tol)
%LPCRA2AR Convert inverse filter autocorrelation coefs to AR filter. AR=(RA)
%
% Usage: (1) ar0=poly([0.5 0.2]);  % ar0=[1 -0.7 0.1]
%            ra0=lpcar2ra(ar0);    % ra0=[1.5 -0.77 0.1]
%            ar1=lpcra2ar(ra0);    % ar1 = ar0
%
%        (2) ar0=poly([0.5 0.2]);                   % ar0=[1 -0.7 0.1]
%            arx=xcorr(ar0);                        % arx=[0.1 -0.77 1.5 -0.77 0.1]
%            ar1=lpcra2ar(arx(length(ar0):end));    % ar1 = ar0
%
%  Inputs: ra(n,p+1) each row is the second half of the autocorrelation of
%                    the coefficients of a stable AR filter of order p
%          tol       tolerance relative to ra(1) [1e-8]
%
% Outputs: ar(n,p+1) each row gives coefficients of an AR filter of order p
%
% This routine uses a Newton-Raphson iteration described in [1] to invert
% the cross-correlation function (as in the second usage example above).
%
% Refs:
% [1]  G. Wilson. Factorization of the covariance generating function of a pure moving average process.
%      SIAM Journal on Numerical Analysis, 6 (1): 1–7, 1969.

%      Copyright (C) Mike Brookes 2015
%      Version: $Id: lpcra2ar.m 6411 2015-07-16 14:38:19Z dmb $
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
persistent p0 i1 i2 j1 j2
if nargin<2
    tol=1e-8;  % tolerance in ra/ra(1)
end
imax=20;
[nf,pp]=size(ra);
if ~numel(p0) || p0~=pp % create index lists for hankel and toeplitz matrices
    p0=pp;
    ix=zeros(1,(pp)*(pp+1)/2);
    nn=1:pp;
    ix(1+(nn-1).*nn/2)=1;
    j1=cumsum(ix);
    i1=cumsum(pp-1+(j1*(1-pp)+pp).*ix)-pp+1;
    j2=pp+1-j1;
    i2=cumsum(((pp+1)*j1-1).*ix-pp-1)+pp+1;
end
ar=zeros(nf,pp);    % space for output filter coefficients
t1=zeros(pp,pp);    % space for hankel coefficient matrix
t2=t1;              % space for toeplitz lower triangular coefficient matrix
ax0=zeros(1,pp);    % temporary filter coefficient row vector
for n=1:nf          % process the input matrix one row at a time
    xa=ra(n,:);     % pick out row n
    ax=ax0;         % initialize ax to have all roots at zero
    ax(1)=sqrt(xa(1)+2*sum(xa(2:end)));
    i=imax;         % maximum number of iterations
    while(i>0)
        t1(i1)=ax(j1); % t1=hankel(ax)
        t2(i2)=ax(j2); % t2=toeplitz(ax,[ax(1) zeros(1,p)])
        ct=ax*t1;
        ax=(xa+ct)/(t1+t2);
        i=min(i-1,i*(max(abs(ct-xa))>tol*xa(1))+1); % do one final iteration after tolerance reached
    end
    ar(n,:)=ax;
end
