function [f,c] = bark2frq(b,m)
%BARK2FRQ  Convert the BARK frequency scale to Hertz FRQ=(BARK)
%
% Inputs: b  matrix of frequencies in Bark
%         m  mode options
%            'h'   use high frequency correction from [1]
%            'l'   use low frequency correction from [1]
%            'H'   do not apply any high frequency correction
%            'L'   do not apply any low frequency correction
%            'u'   unipolar version: do not force b to be an odd function
%                  This has no effect on the default function which is odd anyway
%            's'   use the expression from Schroeder et al. (1979)
%            'g'   plot a graph
%
% Outputs: f  frequency values in Hz
%          c  Critical bandwidth: d(freq)/d(bark)

%   The Bark scale was defined by an ISO committee and published in [2]. It
%   was based on a varienty of experiments on the thresholds for complex
%   sounds, masking, perception of phase and the loudness of complex
%   sounds. The Bark scale is named in honour of Barkhausen, the creator
%   of the unit of loudness level [2]. Critical band k extends
%   from bark2frq(k-1) to bark2frq(k). The inverse function is frq2bark.
%
%   There are many published formulae approximating the Bark scale.
%   The default is the one from [1] but with a correction at high and
%   low frequencies to give a better fit to [2] with a continuous derivative
%   and ensure that 0 Hz = 0 Bark.
%   The h and l mode options apply the corrections from [1] which are
%   not as good and do not give a continuous derivative. The H and L
%   mode options suppress the correction entirely to give a simple formula.
%   The 's' option uses the less accurate formulae from [3] which have been
%   widely used in the lterature.
%
%   [1] H. Traunmuller, Analytical Expressions for the
%       Tonotopic Sensory Scale”, J. Acoust. Soc. Am. 88,
%       1990, pp. 97-100.
%   [2] E. Zwicker, Subdivision of the audible frequency range into
%       critical bands, J Accoust Soc Am 33, 1961, p248.
%   [3] M. R. Schroeder, B. S. Atal, and J. L. Hall. Optimizing digital
%       speech coders by exploiting masking properties of the human ear.
%       J. Acoust Soc Amer, 66 (6): 1647–1652, 1979. doi: 10.1121/1.383662.

%      Copyright (C) Mike Brookes 2006-2010
%      Version: $Id: bark2frq.m 4501 2014-04-24 06:28:21Z dmb $
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
persistent A B C E D P Q R S T U V W X Y Z
if isempty(P)
    A=26.81;
    B=1960;
    C=-0.53;
    E = A+C;
    D=A*B;
    P=(0.53/(3.53)^2);
    V=3-0.5/P;
    W=V^2-9;
    Q=0.25;
    R=20.4;
    xy=2;
    S=0.5*Q/xy;
    T=R+0.5*xy;
    U=T-xy;
    X = T*(1+Q)-Q*R;
    Y = U-0.5/S;
    Z=Y^2-U^2;
end
if nargin<2
    m=' ';
end
if any(m=='u')
    a=b;
else
    a=abs(b);
end
if any(m=='s')
    f=650*sinh(a/7);
else
    if any(m=='l')
        m1=(a<2);
        a(m1)=(a(m1)-0.3)/0.85;
    elseif ~any(m=='L')
        m1=(a<3);
        a(m1)=V+sqrt(W+a(m1)/P);
    end
    if any(m=='h')
        m1=(a>20.1);
        a(m1)=(a(m1)+4.422)/1.22;
    elseif ~any(m=='H')
        m2=(a>X);
        m1=(a>U) & ~m2;
        a(m2)=(a(m2)+Q*R)/(1+Q);
        a(m1)=Y+sqrt(Z+a(m1)/S);
    end
    f=(D*(E-a).^(-1)-B);
end
if ~any(m=='u')
    f=f.*sign(b);          % force to be odd
end
if nargout>1
    [bx,c] = frq2bark(f,m);
end
if ~nargout || any(m=='g')
    [bx,c] = frq2bark(f,m);
    subplot(212)
    semilogy(b,c,'-r');
    ha=gca;
    xlabel('Bark');
    ylabel(['Critical BW (' yticksi 'Hz)']);
    subplot(211)
    plot(b,f,'x-b');
    hb=gca;
    xlabel('Bark');
    ylabel(['Frequency (' yticksi 'Hz)']);
    linkaxes([ha hb],'x');
end