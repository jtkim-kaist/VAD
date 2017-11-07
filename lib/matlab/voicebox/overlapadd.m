function [x,zo]=overlapadd(f,win,inc)
%OVERLAPADD join overlapping frames together X=(F,WIN,INC)
%
% Usage for frequency-domain processing:
%       S=...;                              % input signal
%       OV=2;                               % overlap factor of 2 (4 is also often used)
%       INC=20;                             % set frame increment in samples
%       NW=INC*OV;                          % DFT window length
%       W=sqrt(hamming(NW,'periodic'));     % omit sqrt if OV=4
%       W=W/sqrt(sum(W(1:INC:NW).^2));      % normalize window
%       F=rfft(enframe(S,W,INC),NW,2);      % do STFT: one row per time frame, +ve frequencies only
%       ... process frames ...
%       X=overlapadd(irfft(F,NW,2),W,INC);  % reconstitute the time waveform (omit "X=" to plot waveform)
%
% Inputs:  F(NR,NW) contains the frames to be added together, one
%                   frame per row.
%          WIN(NW)  contains a window function to multiply each frame.
%                   WIN may be omitted to use a default rectangular window
%                   If processing the input in chunks, WIN should be replaced by
%                   ZI on the second and subsequent calls where ZI is the saved
%                   output state from the previous call.
%          INC      gives the time increment (in samples) between
%                   succesive frames [default = NW].
%
% Outputs: X(N,1) is the output signal. The number of output samples is N=NW+(NR-1)*INC.   
%          ZO     Contains the saved state to allow a long signal
%                 to be processed in chunks. In this case X will contain only N=NR*INC
%                 output samples. 
%

%	   Copyright (C) Mike Brookes 2009
%      Version: $Id: overlapadd.m 2470 2012-11-02 15:27:24Z dmb $
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
[nr,nf]=size(f);            % number of frames and frame length
if nargin<2
    win=nf;                 % default increment
end
if isstruct(win)
    w=win.w;
    if ~numel(w) && length(w)~=nf
        error('window length does not match frames size');
    end
    inc=win.inc;
    xx=win.xx;
else
    if nargin<3
        inc=nf;
    end
    if numel(win)==1 && win==fix(win) && nargin<3       % win has been omitted
        inc=win;
        w=[];
    else
        w=win(:).';
        if length(w)~=nf
            error('window length does not match frames size');
        end
        if all(w==1)
            w=[];
        end
    end
    xx=[];      % partial output from previous call is null
end
nb=ceil(nf/inc);        % number of overlap buffers
no=nf+(nr-1)*inc;       % buffer length
z=zeros(no,nb);                      % space for overlapped output speech
if numel(w)
    z(repmat(1:nf,nr,1)+repmat((0:nr-1)'*inc+rem((0:nr-1)',nb)*no,1,nf))=f.*repmat(w,nr,1);
else
    z(repmat(1:nf,nr,1)+repmat((0:nr-1)'*inc+rem((0:nr-1)',nb)*no,1,nf))=f;
end
x=sum(z,2);
if ~isempty(xx)
    x(1:length(xx))=x(1:length(xx))+xx;     % add on leftovers from previous call
end
if nargout>1            % check if we want to preserve the state
    mo=inc*nr;          % completed output samples
    if no<mo
        x(mo,1)=0;
        zo.xx=[];
    else
        zo.xx=x(mo+1:end);
        zo.w=w;
        zo.inc=inc;
        x=x(1:mo);
    end
elseif ~nargout
    if isempty(xx)
        k1=nf-inc;  % dubious samples at start
    else
        k1=0;
    end
    k2=nf-inc;      % dubious samples at end
    plot(1+(0:nr-1)*inc,x(1+(0:nr-1)*inc),'>r',nf+(0:nr-1)*inc,x(nf+(0:nr-1)*inc),'<r', ...
        1:k1+1,x(1:k1+1),':b',k1+1:no-k2,x(k1+1:end-k2),'-b',no-k2:no,x(no-k2:no),':b');
    xlabel('Sample Number');
    title(sprintf('%d frames of %d samples with %.0f%% overlap = %d samples',nr,nf,100*(1-inc/nf),no));
end