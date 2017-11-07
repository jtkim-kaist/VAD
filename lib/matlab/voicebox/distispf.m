function d=distispf(pf1,pf2,mode)
%DISTISPF calculates the Itakura-Saito spectral distance between power spectra D=(PF1,PF2,MODE)
%
% Inputs: PF1,PF2     Power spectra to be compared. Each row represents a power spectrum: the first
%                     and last columns represent the DC and Nyquist terms respectively.
%                     PF1 and PF2 must have the same number of columns.
%
%         MODE        Character string selecting the following options:
%                         'x'  Calculate the full distance matrix from every row of PF1 to every row of PF2
%                         'd'  Calculate only the distance between corresponding rows of PF1 and PF2
%                              The default is 'd' if PF1 and PF2 have the same number of rows otherwise 'x'.
%           
% Output: D           If MODE='d' then D is a column vector with the same number of rows as the shorter of PF1 and PF2.
%                     If MODE='x' then D is a matrix with the same number of rows as PF1 and the same number of columns as PF2'.
%
% The Itakura-Saito spectral distance is the average over +ve and -ve frequency of 
%
%                      pf1/pf2 - log(pf1/pf2) - 1     =     exp(v) - v - 1         where v=log(pf1/pf2)
%
% The Itakura-Saito distance is asymmetric: pf1>pf2 contributes more to the distance than pf2>pf1. 
% A symmetrical version is the COSH distance: distchpf(x,y)=(distispf(x,y)+distispf(y,x))/2

% The Itakura-Saito distance can also be calculated directly from AR coefficients; providing np is large
% enough, the values of d0 and d1 in the following will be very similar:
%
%         np=255; d0=distisar(ar1,ar2); d1=distispf(lpcar2pf(ar1,np),lpcar2pf(ar2,np))
%

% Ref: A.H.Gray Jr and J.D.Markel, "Distance measures for speech processing", IEEE ASSP-24(5): 380-391, Oct 1976
%      L. Rabiner abd B-H Juang, "Fundamentals of Speech Recognition", Section 4.5, Prentice-Hall 1993, ISBN 0-13-015157-2
%      F.Itakura & S.Saito, "A statistical method for estimation of speech spectral density and formant frequencies",
%                            Electronics & Communications in Japan, 53A: 36-43, 1970.


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: distispf.m 713 2011-10-16 14:45:43Z dmb $
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

[nf1,p2]=size(pf1);
p1=p2-1;
nf2=size(pf2,1);
if nargin<3 | isempty(mode) mode='0'; end
if any(mode=='d') | (mode~='x' & nf1==nf2)
   nx=min(nf1,nf2);
   r=pf1(1:nx,:)./pf2(1:nx,:);
   q=r-log(r);
   d=(sum(q(:,2:p1),2)+0.5*(q(:,1)+q(:,p2)))/p1-1;
else
   r=permute(pf1(:,:,ones(1,nf2)),[1 3 2])./permute(pf2(:,:,ones(1,nf1)),[3 1 2]);
   q=r-log(r);
   d=(sum(q(:,:,2:p1),3)+0.5*(q(:,:,1)+q(:,:,p2)))/p1-1;
end
