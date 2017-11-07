function f=midi2frq(n,s)
%MIDI2FRQ	Convert musical note numbers to frequencies F=(N,S)
%		s is:	'e' equal tempered (default)
%			'p' pythagorean scale
%			'j' just intonation
%
% notes are numbered in semitones with middle C being 60
% On the equal tempered scale, note 69 (the A above middle C)
% has a frequency of 440 Hz.
%
% see FRQ2NOTE for the inverse transform

% Pythagorean
%     sharps 1 2187/2048 9/8 19683/16384 81/64 4/3 729/512  3/2 6561/4096 27/16 59049/32768 243/128 2
%     flats  1 256/243   9/8 32/27       81/64 4/3 1024/729 3/2 128/81    27/16 16/9        243/128 2
%
% Just Intonation
%     sharps 1 25/24 9/8 75/64 5/4 4/3 45/32  3/2 25/16 5/3 225/128 15/8 2
%     flats  1 16/15 9/8 6/5   5/4 4/3 108/75 3/2 8/5   5/3 18/10   15/8 2


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: midi2frq.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin > 1
  if s(1)=='p'
    r=[256/243 9/8 32/27 81/64 4/3 729/512 3/2 128/81 27/16 16/9 243/128];
  elseif s(1)=='j'
    r=[16/15 9/8 6/5 5/4 4/3 36/25 3/2 8/5 5/3 9/5 15/8];
  else
    r=0;
  end
  if r(1)
    c=[0 0 12*log(r)/log(2)-(1:11) 0];
    nm=mod(n,12);
    na=floor(nm);
    nb=nm-na;
    f=440*exp((n+c(na+2).*(1-nb)+c(na+3).*nb-69)*log(2)/12);
  else
    f=440*exp((n-69)*log(2)/12);
  end
else
  f=440*exp((n-69)*log(2)/12);
end

