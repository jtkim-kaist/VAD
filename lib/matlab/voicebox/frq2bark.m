function [b,c] = frq2bark(f,m)
%FRQ2BARK  Convert Hertz to BARK frequency scale BARK=(FRQ)
%       bark = frq2bark(frq) converts a vector of frequencies (in Hz)
%       to the corresponding values on the BARK scale.
% Inputs: f  matrix of frequencies in Hz
%         m  mode options
%            'h'   use high frequency correction from [1]
%            'l'   use low frequency correction from [1]
%            'H'   do not apply any high frequency correction
%            'L'   do not apply any low frequency correction
%            'z'   use the expressions from Zwicker et al. (1980) for b and c
%            's'   use the expression from Schroeder et al. (1979)
%            'u'   unipolar version: do not force b to be an odd function
%                  This has no effect on the default function which is odd anyway
%            'g'   plot a graph
%
% Outputs: b  bark values
%          c  Critical bandwidth: d(freq)/d(bark)

%   The Bark scale was defined by in ISO532 and published in [2]. It
%   was based on a varienty of experiments on the thresholds for complex
%   sounds, masking, perception of phase and the loudness of complex
%   sounds. The Bark scale is named in honour of Barkhausen, the creator
%   of the unit of loudness level [2]. Frequency f lies in critical 
%   band ceil(frq2bark(f)). The inverse function is bark2frq.
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
%   The 'z' option uses the formulae from [4] in which the c output
%   is not exactly the reciprocal of the derivative of the bark function.
%
%   [1] H. Traunmuller, Analytical Expressions for the
%       Tonotopic Sensory Scale”, J. Acoust. Soc. Am. 88,
%       1990, pp. 97-100.
%   [2] E. Zwicker, Subdivision of the audible frequency range into
%       critical bands, J Accoust Soc Am 33, 1961, p248.
%   [3] M. R. Schroeder, B. S. Atal, and J. L. Hall. Optimizing digital
%       speech coders by exploiting masking properties of the human ear.
%       J. Acoust Soc Amer, 66 (6): 1647–1652, 1979. doi: 10.1121/1.383662.
%   [4] E. Zwicker and E. Terhardt.  Analytical expressions for
%       critical-band rate and critical bandwidth as a function of frequency.
%       J. Acoust Soc Amer, 68 (5): 1523–1525, Nov. 1980.

%   The following code reproduces the graphs 3(c) and 3(d) from [1].
%       b0=(0:0.5:24)';
%       f0=[[2 5 10 15 20 25 30 35 40 45 51 57 63 70 77 ...
%           84 92 100 108 117 127 137 148 160 172 185 200 ...
%           215 232 250 270 290 315]*10 [34 37 40 44 48 53 ...
%           58 64 70 77 85 95 105 120 135 155]*100]';
%       b1=frq2bark(f0);      b2=frq2bark(f0,'lh');
%       b3=frq2bark(f0,'LH'); b4=frq2bark(f0,'z');
%       plot(b0,[b0 b1 b2 b3 b4]-repmat(b0,1,5));
%       xlabel('Frequency (Bark)'); ylabel('Error (Bark)');
%       legend('Exact','voicebox','Traunmuller1990', ...
%              'Traunmuller1983','Zwicker1980','Location','South');

%      Copyright (C) Mike Brookes 2006-2010
%      Version: $Id: frq2bark.m 4501 2014-04-24 06:28:21Z dmb $
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
persistent A B C D P Q R S T U
if isempty(P)
    A=26.81;
    B=1960;
    C=-0.53;
    D=A*B;
    P=(0.53/(3.53)^2);
    Q=0.25;
    R=20.4;
    xy=2;
    S=0.5*Q/xy;
    T=R+0.5*xy;
    U=T-xy;
end
if nargin<2
    m=' ';
end
if any(m=='u')
    g=f;
else
    g=abs(f);
end
if any(m=='z')
    b=13*atan(0.00076*g)+3.5*atan((f/7500).^2);
    c=25+75*(1+1.4e-6*f.^2).^0.69;
elseif any(m=='s')
    b=7*log(g/650+sqrt(1+(g/650).^2));
    c=cosh(b/7)*650/7;
else
    b=A*g./(B+g)+C;
    d=D*(B+g).^(-2);
    if any(m=='l')
        m1=(b<2);
        d(m1)=d(m1)*0.85;
        b(m1)=0.3+0.85*b(m1);
    elseif ~any(m=='L')
        m1=(b<3);
        b(m1)=b(m1)+P*(3-b(m1)).^2;
        d(m1)=d(m1).*(1-2*P*(3-b(m1)));
    end
    if any(m=='h')
        m1=(b>20.1);
        d(m1)=d(m1)*1.22;
        b(m1)=1.22*b(m1)-4.422;
    elseif ~any(m=='H')
        m2=(b>T);
        m1=(b>U) & ~m2;
        b(m1)=b(m1)+S*(b(m1)-U).^2;
        b(m2)=(1+Q)*b(m2)-Q*R;
        d(m2)=d(m2).*(1+Q);
        d(m1)=d(m1).*(1+2*S*(b(m1)-U));
    end
    c=d.^(-1);
end
if ~any(m=='u')
    b=b.*sign(f);          % force to be odd
end

if ~nargout || any(m=='g')
    subplot(212)
    semilogy(f,c,'-r');
    ha=gca;
    ylabel(['Critical BW (' yticksi 'Hz)']);
    xlabel(['Frequency (' xticksi 'Hz)']);
    subplot(211)
    plot(f,b,'x-b');
    hb=gca;
    ylabel('Bark');
    xlabel(['Frequency (' xticksi 'Hz)']);
    linkaxes([ha hb],'x');
end
