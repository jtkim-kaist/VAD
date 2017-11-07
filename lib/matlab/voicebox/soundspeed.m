function [v,d,z]=soundspeed(t,p,m,g)
%SOUNDSPEED gives the speed of sound, density of air and acoustic impedance as a function of temp & pressure [V,D,Z]=(T,P,M,G)
%
%  Inputs:  T        air temperature in Celsius  [20 deg C]
%           P        air pressure [1 atm]
%           M        average molecular weight of air [0.0289644 kg/mol]
%           G        adiabatic constant for air [1.4]
%
% Outputs:  V        is the speed of sound in m/s
%           D        density of air in kg/m^3
%           Z        characteristic impedance of air Pa.s/m

% Notes: (1) Sound pressure is often measured in dB (SPL) relative to 20uPa [20*log10(p/p0)]
%            Sound pressure is inversely proportional to distance.
%        (2) Sound intensity is often measured indB relative to pW/m^2,[10*log10(J*10^12)]
%            Intensity is inversely proportional to distance squared.
%        (3) Intensity * impedance = pressure^2, so with the default values, 1 pW/m^2 = sqrt(Z) = 20.33 uPa
%            So: X dB (SPL) = X-93.98 dB (Pa) = X-0.14 dB (pW/m^2) =  X-120.14 dB (W/m^2)
%            where 93.98=20*log10(50000), 0.14=10*log10(Z/400), 120.14=10*log10(1e12*z/400)
%        (4) The default air pressure (which does not affect sound speed) in various units is:
%            1 atm = 101325 Pa = 1.01325 bar = 1.0332 at = 760 torr = 14.696 psi

%	   Copyright (C) Mike Brookes 2006
%      Version: $Id: soundspeed.m 4967 2014-08-05 18:21:35Z dmb $
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

if nargin<4
    g=1.4;
    if nargin<3
        m=0.0289644;      % gm/mol
        if nargin<2
            p=1;
            if nargin<1
                t=20;
            end
        end
    end
end
p=p*101325; % convert pressure: atm to pascal
k=t+273.15;  % absolute temperature
r=8.3144;  % J/(mol K) universal gas constant
d=p*m/(r*k);
v=sqrt(g*r*k/m);
z=v*d;
