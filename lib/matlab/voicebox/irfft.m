function x=irfft(y,n,d)
%IRFFT    Inverse fft of a conjugate symmetric spectrum X=(Y,N,D)
%
% Inputs:  Y(M)   The first half of a complex spectrum
%          N      The number of output points to generate (default: 2M-2)
%          D      The dimension along which to perorm the transform
%                 (default: first non-singleton dimension of Y)
%
% Outputs: X(N)   Real inverse dft of Y
%
% This routine calculates the inverse DFT of a conjugate-symmetric to give a real-valued
% output of dimension N. Only the first half of the spectrum need be supplied: if N is even,
% this includes the Nyquist term and is of dimension M=N/2 + 1 whereas if N is odd then there is
% no Nyquist term and the input is of dimension M=(N+1)/2.
% Note that the default value of N is always even so that N must be given explicitly
% if it is odd.
%
% See also the forward transform: RFFT

%      Copyright (C) Mike Brookes 2009
%      Version: $Id: irfft.m 2460 2012-10-29 22:20:45Z dmb $
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

s=size(y);
ps=prod(s);
ns=length(s);
if ps==1
    x=y;
else
    if nargin <3 || isempty(d)
        d=find(s>1,1);
    end
    m=s(d);
    k=ps/m;     % number of fft's to do
    if d==1
        v=reshape(y,m,k);
    else
        v=reshape(permute(y,[d:ns 1:d-1]),m,k);
    end
    if nargin<2 || isempty(n)
        n=2*m-2;        % default output length
    else
        mm=1+fix(n/2);          % expected input length
        if mm>m v=[v; zeros(mm-m,k)];   % zero pad
        elseif mm<m v(mm+1:m,:)=[];     % or truncate
        end
        m=mm;
    end
    if rem(n,2)		% odd output length
        x=real(ifft([v;conj(v(m:-1:2,:))],[],1));    % do it the long way
    else			% even output length
        v(m,:)=real(v(m,:));	% force nyquist element real
        w=ones(1,k);
        %  t=[cumprod([-0.5i; exp(2i*pi/n)*ones(m-2,1)]); 0.5i];
        t=-0.5i* exp((2i*pi/n)*(0:m-1)).';
        z=(t(:,w)+0.5).*(conj(flipud(v))-v)+v;
        z(m,:)=[];
        zz=ifft(z,[],1);
        x=zeros(n,k);
        x(1:2:n,:)=real(zz);
        x(2:2:n,:)=imag(zz);
    end
    s(d)=n;         % change output dimension
    if d==1
        x=reshape(x,s);
    else
        x=permute(reshape(x,s([d:ns 1:d-1])),[ns+2-d:ns 1:ns+1-d]);
    end
end
