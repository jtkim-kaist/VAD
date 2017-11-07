function [n,t]=frq2midi(f)
%FRQ2MIDI Convert frequencies to musical note numbers [N,T]=(F)
% notes are numbered in semitones with middle C being 60
% Note 69 (the A above middle C) has a frequency of 440 Hz.
% These note numbers are used by MIDI. Note numbers are not necessarily
% integers.
%
% t is a text representation of the note in which
% C4# denotes C sharp in octave 4. Octave 4 goes
% from middle C up to the B above middle C. For the white
% notes on the piano, the third character is a space.
%
% Negative frequencies are equivalent to positive frequencies
% except that flats will be used instead of sharps. Thus
% C4# would become D4-
%
% see MIDI2FRQ for the inverse transform




%      Copyright (C) Mike Brookes 1998
%      Version: $Id: frq2midi.m 713 2011-10-16 14:45:43Z dmb $
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

n=(69+12*log(abs(f)/440)/log(2));
if nargout > 1
  m=round(n(:));
  o=floor(m/12)-1;
  m=m-12*o+6*sign(f(:))-5;
  a=('CDDEEFGGAABBCCDDEFFGGAAB')';
  b=(' - -  - - -  # #  # # # ')';
  t=setstr([a(m) mod(o,10)+'0' b(m)]);
end
