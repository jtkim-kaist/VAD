function d=distisar(ar1,ar2,mode)
%DISTISAR calculates the Itakura-Saito distance between AR coefficients D=(AR1,AR2,MODE)
%
% Inputs: AR1,AR2     AR coefficient sets to be compared. Each row contains a set of coefficients.
%                     AR1 and AR2 must have the same number of columns.
%
%         MODE        Character string selecting the following options:
%                         'x'  Calculate the full distance matrix from every row of AR1 to every row of AR2
%                         'd'  Calculate only the distance between corresponding rows of AR1 and AR2
%                              The default is 'd' if AR1 and AR2 have the same number of rows otherwise 'x'.
%           
% Output: D           If MODE='d' then D is a column vector with the same number of rows as the shorter of AR1 and AR2.
%                     If MODE='x' then D is a matrix with the same number of rows as AR1 and the same number of columns as AR2'.
%
% The Itakura-Saito spectral distance is the average over +ve and -ve frequency of 
%
%                      pf1/pf2 - log(pf1/pf2) - 1     =     exp(v) - v - 1         where v=log(pf1/pf2)
%
% The Itakura-Saito distance is asymmetric: pf1>pf2 contributes more to the distance than pf2>pf1. 
% A symmetrical version is the COSH distance: distchpf(x,y)=(distispf(x,y)+distispf(y,x))/2
%
% The I-S distance can be expressed as ar2*toeplitz(lpcar2rr(ar1))*ar2' + log((ar1(1)/ar2(1)).^2) - 1
% but this is not how we actually calculate it.


% Since the power spectrum is the fourier transform of the autocorrelation, we can calculate
% the average value of p1/p2 by taking the 0'th order term of the convolution of the autocorrelation
% functions associated with p1 and 1/p2. Since 1/p2 corresponds to an FIR filter, this convolution is
% a finite sum even though the autocorrelation function of p1 is infinite in extent.
% The average value of log(pf1) is equal to log(ar1(1)^-2) where ar1(1) is the 0'th order AR coefficient.

% The Itakura-Saito distance can also be calculated directly from the power spectra; providing np is large
% enough, the values of d0 and d1 in the following will be very similar:
%
%         np=255; d0=distisar(ar1,ar2); d1=distispf(lpcar2pf(ar1,np),lpcar2pf(ar2,np))
%
% Autocorrelation LPC analysis is equivalent to minimizing the Itakura-Saito difference between the
% signal spectrum and that of the all-pole LPC filter, i.e. distispf(pf,lpcar2pf(ar0,np)).
% Moreover, if ar0 is the LPC filter and ar is any  other all-pole filter, the I-S distance has the
% following additive property:
%
%               distispf(pf,lpcar2pf(ar,np)) = distispf(pf,lpcar2pf(ar0,np)) + distisar(ar0,ar)

% Ref: A.H.Gray Jr and J.D.Markel, "Distance measures for speech processing", IEEE ASSP-24(5): 380-391, Oct 1976
%      L. Rabiner abd B-H Juang, "Fundamentals of Speech Recognition", Section 4.5, Prentice-Hall 1993, ISBN 0-13-015157-2
%      F.Itakura & S.Saito, "A statistical method for estimation of speech spectral density and formant frequencies",
%                            Electronics & Communications in Japan, 53A: 36-43, 1970.

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: distisar.m 713 2011-10-16 14:45:43Z dmb $
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

[nf1,p1]=size(ar1);
nf2=size(ar2,1);
m2=lpcar2ra(ar2);
m2(:,1)=m2(:,1)*0.5;
if nargin<3 | isempty(mode) mode='0'; end
if any(mode=='d') | (mode~='x' & nf1==nf2)
   nx=min(nf1,nf2);
   d=2*sum(lpcar2rr(ar1(1:nx,:)).*m2(1:nx,:),2)-log((ar2(1:nx,1)./ar1(1:nx,1)).^2)-1;;
else
   d=2*lpcar2rr(ar1)*m2'-log((ar1(:,1).^(-1)*ar2(:,1)').^2)-1;
end
