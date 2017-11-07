function x=lambda2rgb(l,m)
%LAMBDA2XYZ Convert wavelength to XYZ or RGB colour space X=(L,M)
%
%Inputs:   l(n,1)  column vector of wavelengths in nanometres
%          m       mode:
%                    r - output is [R G B] using the 1931 observer with negative values eliminated [default]
%                    s - output is [R G B] using the 1931 observer with signed values
%                    x - output ix [X Y Z] using the 1931 observer
%
%                    Use uppercase 'X', 'R' etc for the 1964 observer instead
%
% Outputs: x(n,3)  tristimulus output values
%
% The formulae are taken from [1] and were obtained by numerical curve
% fitting to the CIE standard observer data available from [2].
%
% References:
%  [1]	C. Wyman, P.-P. Sloan, and P. Shirley.
%       Simple analytic approximations to the CIE XYZ color matching functions.
%       Journal of Computer Graphics Techniques, 2 (2): 1–11, 2013.
%  [2]	D. Wyble. Useful color data. Website, Rochester Institute of Technology, 2001.
%       URL http://www.rit.edu/cos/colorscience/rc_useful_data.php

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: lambda2rgb.m 4899 2014-07-23 08:31:57Z dmb $
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
persistent c d xr rx
if isempty(c)
    c=[1.065 -0.5/33.33^2 595.8 0.366 -0.5/19.44^2 446.8 1.014 -0.5/0.075^2 log(556.3) 1.839 -0.5/0.051^2 log(449.8)];
    d=[0.398 -1250 -570.1 -log(1014) 1.132 -234 1338 1300 -log(743.5) 1.011 -0.5/46.14^2 556.1 2.06 -32 265.8 -log(180.4)];
    xr=[0.49 0.31 0.2; 0.17697 0.8124 0.01063; 0 0.01 0.99];
    xr=xr'/xr(2);
    rx=inv(xr);
end
if nargin<1
    l=[];
end
lv=l(:);
if nargin<2
    m='r';
end
lm=lower(m);
if m==lm % use 1931 standard observer
    ll=log(lv);
    x=[c(1)*exp(c(2)*(lv-c(3)).^2)+c(4)*exp(c(5)*(lv-c(6)).^2) c(7)*exp(c(8)*(ll-c(9)).^2) c(10)*exp(c(11)*(ll-c(12)).^2)];
else % use 1964 standard observer
    x=[d(1)*exp(d(2)*(log(lv-d(3))+d(4)).^2)+d(5)*exp(d(6)*(log(d(7)-min(lv,d(8)))+d(9)).^2) d(10)*exp(d(11)*(lv-d(12)).^2) d(13)*exp(d(14)*(log(lv-d(15))+d(16)).^2)];
end
if lm=='s'
    x=x*rx;
elseif lm=='r'
    x=max(x*rx,0);
end
if ~nargout
    if numel(l)<10
        la=linspace(360,740,200)';
        xa=lambda2rgb(la,m);
    else
        la=l;
        xa=x;
    end
    plot(la,xa(:,1),'r-',la,xa(:,2),'g-',la,xa(:,3),'b-');
    if numel(l)<10 && numel(l)>0
        hold on
        plot(l,x(:,1),'ro',l,x(:,2),'go',l,x(:,3),'bo');
        hold off
    end
    xlabel('Wavelength, \lambda (nm)');
    switch lm
        case 'x'
            cstr='XYZ';
        case 's'
            cstr='RGB';
        case 'r'
            cstr='RGB+';
    end
    if lm==m
        yr=1931;
    else
        yr=1964;
    end
    legend(cstr(1), cstr(2), cstr(3));
    title(sprintf('%d Standard Observer (%s)',yr,cstr));
end


