function u=glotlf(d,t,p)
%GLOTLF   Liljencrants-Fant glottal model U=(D,T,P)
% d is derivative of flow waveform: must be 0, 1 or 2
% t is a vector of time instants at which to calculate the
%   waveform. Time units are in fractions of a cycle.
% p is a vector of 3 parameters defining the waveform
%    p(1) is the time at which ugd has its peak negative value. This we define as the
%         start of the closed phase. p(1) is therefore the open/closed interval ratio.
%    p(2) is the reciprocal of the peak negative value of ugd(t)
%    p(3) is the fraction of the open phase for which ugd(t) is negative. That is, it is
%         it is the time between the peak flow and the end of the open phase expressed
%         as a fraction of the open phase.
%
% Note: this signal has not been low-pass filtered
% and will therefore be aliased [this is a bug]
%
% Usage example:	ncyc=5;
%			period=80;
%			t=0:1/period:ncyc;
%			ug=glotlf(0,t);
%			plot(t,ug)


%      Copyright (C) Mike Brookes 1998
%      Version: $Id: glotlf.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin < 2
  tt=(0:99)'/100;
else
  tt=t-floor(t);
end
u=zeros(size(tt));
de=[0.6 0.1 0.2]';
if nargin < 3
  p=de;
elseif length(p)<2
  p=[p(:); de(length(p)+1:2)];
end

% Calculate the parameters in terms of ugd(t), the glottal flow derivative

te=p(1);            % ugd(te) is the negative peak.
mtc=te-1;
e0=1;
wa=pi/(te*(1-p(3)));
a=-log(-p(2)*sin(wa*te))/te;
inta=e0*((wa/tan(wa*te)-a)/p(2)+wa)/(a^2+wa^2);

% if inta<0 we should reduce p(2)
% if inta>0.5*p(2)*(1-te) we should increase p(2)

rb0=p(2)*inta;
rb=rb0;

% Use Newton to determine closure time constant
% so that flow starts and ends at zero.

for i=1:4
  kk=1-exp(mtc/rb);
  err=rb+mtc*(1/kk-1)-rb0;
  derr=1-(1-kk)*(mtc/rb/kk)^2;
  rb=rb-err/derr;
end
e1=1/(p(2)*(1-exp(mtc/rb)));


ta=tt<te;
tb=~ta;

if d==0
  u(ta)=e0*(exp(a*tt(ta)).*(a*sin(wa*tt(ta))-wa*cos(wa*tt(ta)))+wa)/(a^2+wa^2);
  u(tb)=e1*(exp(mtc/rb)*(tt(tb)-1-rb)+exp((te-tt(tb))/rb)*rb);
elseif d==1
  u(ta)=e0*exp(a*tt(ta)).*sin(wa*tt(ta));
  u(tb)=e1*(exp(mtc/rb)-exp((te-tt(tb))/rb));
elseif d==2
  u(ta)=e0*exp(a*tt(ta)).*(a*sin(wa*tt(ta))+wa*cos(wa*tt(ta)));
  u(tb)=e1*exp((te-tt(tb))/rb)/rb;
else
  error('Derivative must be 0,1 or 2');
end
