function [b,a]=potsband(fs)
%POTSBAND Design filter for 300-3400 telephone bandwidth [B,A]=(FS)
%
%Input: FS=sample frequency in Hz
%
%Output: B/A is a discrete time bandpass filter with a passband gain of 1
%
%The filter meets the specifications of G.151 for any sample frequency
%and has a gain of -3dB at the passband edges.


%      Copyright (C) Mike Brookes 1998
%      Version: $Id: potsband.m 713 2011-10-16 14:45:43Z dmb $
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

szp=[0.19892796195357i; -0.48623571568937+0.86535995266875i]; 
szp=[[0; -0.97247143137874] szp conj(szp)];
% s-plane zeros and poles of high pass 3'rd order chebychev2 filter with -3dB at w=1
zl=2./(1-szp*tan(300*pi/fs))-1;
al=real(poly(zl(2,:)));
bl=real(poly(zl(1,:)));
sw=[1;-1;1;-1];
bl=bl*(al*sw)/(bl*sw);
zh=2./(szp/tan(3400*pi/fs)-1)+1;
ah=real(poly(zh(2,:)));
bh=real(poly(zh(1,:)));
bh=bh*sum(ah)/sum(bh);
b=conv(bh,bl);
a=conv(ah,al);