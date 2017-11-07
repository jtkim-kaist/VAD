function [u,v,w]=sphrharm(m,a,b,c,d)
%SPHRHARM  forward and inverse spherical harmonic transform
%
% Usage: (1) y=('f',n,x)      % Calculate complex transform of spatial data x up to order n
%
%        (2) y=('fr',n,x)     % Calculate real transform of spatial data x(ne,na) up to order n
%                             % x is specified on a uniform grid with inclination e=(0.5:ne-0.5)*pi/ne
%                             % and azimuth a=(0:na-1)*2*pi/na. The North pole has inclination 0 and
%                             % azimuths increase going East.
%
%        (3) y=('fd',n,e,a,v) % Calculate transform of impulse at (e(i),a(i)) of height v(i) up to order n
%                             % e(i), a(i) and v(i) should be the same dimension; v defaults to 1.
%
%        (4) x=('i',y,ne,na)  % Calculate spatial data from spherical harmonics y on
%                             % a uniform grid with ne inclinations and na azimuths
%
%        (5) [e,a,w]=('cg',ne,na)   % Calculate the inclinations, azimuths and
%                                   % quadrature weights of a Gaussian sampling grid
%
%        (6) n=2;             % illustrate real basis functions upto order 2
%            for i=0:n
%              for j=-i:i
%                subplot(n+1,2*n+1,i*(2*n+1)+j+n+1);
%                sphrharm('irp',[zeros(1,i^2+i+j) 1],25,25);
%              end
%            end
%
%  Inputs:  m string specifying various options:
%               'f' forward transform
%               'i' inverse transform
%               'c' calculate the coordinates of the sampling grid [default if f or i omitted]
%               'r' real spherical harmonics instead of complex
%               'u' uniform inclination grid:  e=(0.5:ne-0.5)*pi/ne [default]
%                      for invertibility, ne>=2N+1 for order N
%               'U' uniform inclination grid:  e=(0:ne-1)*pi/ne
%                      for invertibility, ne>=2N+1 for order N
%               'g' gaussian inclination grid (non-uniform but fewer samples needed)
%                      for invertibility, ne>=N+1 for order N
%               'a' arbitrary (specified by user - inverse transform only)
%               'd' delta function [forward transform only]
%               'p' plot result
%               'P' plot with colourbar
%
%           The remaining inputs depend on the mode specified:
%               'f'  a         order of transform
%                    b(ne,na)  input data array on the chosen inclination-azimuth sampling grid
%               'fd' a         order of transform
%                    b(k)      inclinations of k delta functions
%                    c(k)      azimuths of k delta functions
%                    d(k)      amplitudes of k delta functions [default=1]
%               'i'  a()       spherical harmonics as a single row vector in the order:
%                                 (0,0),(1,-1),(1,0),(1,1),(2,-2),(2,-1),...
%                              To access as a 2-D harmonic: Y(n,m)=a(n^2+n+m+1) where n=0:N, m=-n:n
%                    b         number of inclination values in output grid
%                    c         number of azimuth values in output grid
%               'c'  a         number of inclination values in output grid
%                    b         number of azimuth values in output grid
%
% Outputs:  u output data according to transform direction:
%                'f': spherical harmonics as a single row vector in the order:
%                        (0,0),(1,-1),(1,0),(1,1),(2,-2),(2,-1),...
%                     To access as a 2-D harmonic: Y(n,m)=u(n^2+n+m+1) where n=0:N, m=-n:n
%                'i': u(ne,na) is spatial data sampled on an azimuth-inclination grid
%                     with ne inclination points (in 0:pi) and na azimuth points (in 0:2*pi)
%                'c': u(ne) gives inclination grid positions with 0 being the North pole
%           v(na) gives azimuth grid positions
%           w(ne) gives the quadrature weights used in the transform
%
% Suppose f(e,a) is a complex-valued function defined on the surface of the
% sphere (0<=e<=pi, 0<=a<2pi) where e=inclination=pi/2-elevation and
% a=azimuth. (0,*) is the North pole, (pi/2,0) is on the equator near
% Ghana and (pi/2,pi/2) is on the equator near India.
%
% We can approximate f(e,a) using complex spherical harmonics as
% f(e,a) = sum(F(n,m)*u(n,m,e,a)) where the sum
% is over n=0:N,m=-n:n giving a total of (N+1)^2 coefficients F(n,m).
%
% If f(e,a) happens to be real-valued, then we can transform the
% coefficients using G(n,0)=F(n,0) and, for m>0,
% [G(n,m); G(n,-m)]=sqrt(0.5)*[1 1; 1i -1i]*[F(n,m); F(n,-m)]
% to give the real spherical harmonic coefficients G(n,m).
%
% The basis functions u(n,m,e,a) are orthonormal under the inner product
% <u(e,a),v(e,a)> = integral(u(e,a)*v'(e,a)*sin(e)*de*da).
%
%  Minimum spatial grid for invertibility:
%     Gaussian grid: ne >= N+1,  na >= 2N+1
%     Uniform grid:  ne >= 2N+1, na >= 2N+1
%     An inverse transform followed by a forward transform will restore the
%     original coefficient values provided that the sampling grid is large
%     enough. The advantage of the Gaussian grid is that it can be smaller.
%
%   Data formats:
%     (1) Spatial Data: x=sphrharm('i',y,ne,na)
%            x(1:ne,1:na) is spatial data sampled on an azimuth-inclination grid
%                         with ne inclination points (in 0:pi) and
%                         na azimuth points (in 0:2*pi)
%     (2) Spherical harmonics: y=sphrharm('f',n,x)
%            y(1:(n+1)^2)  spherical harmonics as a single row vector
%                          in the order: (0,0),(1,-1),(1,0),(1,1),(2,-2),(2,-1),...
%                          To access as a 2-D harmonic, use
%                          Y(n,m)=y(n^2+n+m+1)   where n=0:N, m=-i:i
%     (3) Sampling Grid:  [e,a,w]=sphrharm('c',ne,na)
%            e(1:ne)   monotonically increasing inclination values (0<=e<=pi) where 0 is the North pole
%            a(1:na)   monotonically increasing azimuth values (0<=a<2pi)
%            w(1:ne)   Quadrature weights
%

% future options
%    direction: [m=rotation matrix]
%    transform: [h=complex harmonics but only the positive ones specified]
%    grid       d=delta function [forward transform only]
%               [s=sparse azimuth]
%               [z=include both north and south pole]
% 
% bugs:
%   (1) we could save space and time by taking advantage of symmetry of legendre polynomials
%   (2) the number of points in the inclination (elevation) grid must, for now, be even if the 'U' option is chosen
%   (3) should ensure that the negative nyquist azimuth frequency is not used
%   (4) save time by manipulating only the necessary 2*m columns of the da matrix
%   (5) should make non-existant coefficients black in plots
%   (6) using 'surf' for plots adds an offset in azimuth and elevation
%       because colours correspond to vertices not faces
%   (7) the normalization for mode 'fd' seems incorrect
%   (8) mode 'fd' should allow multiple impulses to be summed

% errors:

% tests:
% check for inverse transform n=4; m=4; [ve,va]=spvals(8,8,n,m); ve*va-sphrharm('iur',[zeros(1,n^2+n+m) 1],8,8)
% check inverse followed by forward: no=4; h=rand(1,(no+1)^2); max(abs(sphrharm('fur',sphrharm('iur',h,10,10),no)-h))
% same but complex: no=4; h=rand(1,(no+1)^2); max(abs(sphrharm('fu',sphrharm('iu',h,10,10),no)-h))
% same but gaussian grid: no=4; h=rand(1,(no+1)^2); max(abs(sphrharm('fg',sphrharm('ig',h,6,10),no)-h))

%      Copyright (C) Mike Brookes 2009
%      Version: $Id: sphrharm.m 2467 2012-11-02 08:12:41Z dmb $
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

% decode option string
mv='c u ';       % mv(1)=transform, mv(2)=real/complex, mv(3)=grid, mv(4)=plot
if ~nargout
    mv(4)='P';
end
mc='ficruUgadpP';
mi=[1 1 1 2 3 3 3 3 3 4 4];
for i=1:length(m)
    j=find(m(i)==mc);
    if ~isempty(j)
        mv(mi(j))=m(i);
    end
end

switch mv(1)
    case 'f' % forward transform [sp]=('f',order,data)
        if mv(3)~='d'
            % input data has elevations down the columns and azimuths along the rows
            % output data is in the form:
            [ne,na]=size(b);
            da=fft(b,[],2)*sqrt(pi)/na;
            if mv(2)=='r' % only actually need to do this for min(a,floor(na/2)) values (but nyquist is tricky)
                ix=2:ceil(na/2);
                iy=na:-1:na+2-ix(end);
                da(:,ix)=da(:,ix)+da(:,iy);
                da(:,iy)=1i*(da(:,ix)-2*da(:,iy));  % note
                da(:,1)=da(:,1)*sqrt(2);
            else
                da=da*sqrt(2);
            end
            [ue,we,lgp]=sphrharp(mv(3),ne,a);
            da=da.*repmat(we',1,na);
            u=zeros(1,(a+1)^2);
            i=0;
            for nn=0:a
                % we could vectorize this to avoid the inner loop
                % we should ensure the the negative nyquist value is not used
                for mm=-nn:nn
                    i=i+1;
                    u(i)=lgp(1+nn*(nn+1)/2+abs(mm),:)*da(:,1+mod(mm,na));
                end
            end
        else       % forward transform of impulses [sp]=('fd',order,e,a,value)
            if nargin<5
                d=ones(size(b));
            end
            [ue,we,lgp]=sphrharp('a',b(:),a);
            uu=zeros(1,(a+1)^2);     % reserve space for spherical harmonic coefficients
            u=uu;
            for j=1:numel(b)
                i=0;
                if mv(2)=='r'
                    for nn=0:a
                        % we could vectorize this to avoid the inner loop
                        % we should ensure the the negative nyquist value is not used
                        i=i+1;
                        uu(i+nn)=lgp(1+nn*(nn+1)/2,1);
                        for mm=1:nn
                            uu(i+nn-mm)=lgp(1+nn*(nn+1)/2+abs(mm),j)*sin(mm*c(j))*sqrt(2);
                            uu(i+nn+mm)=lgp(1+nn*(nn+1)/2+abs(mm),j)*cos(mm*c(j))*sqrt(2);
                        end
                        i=i+2*nn;
                    end
                else
                    for nn=0:a
                        % we could vectorize this to avoid the inner loop
                        % we should ensure the the negative nyquist value is not used
                        for mm=-nn:nn
                            i=i+1;
                            uu(i)=lgp(1+nn*(nn+1)/2+abs(mm),j)*exp(-1i*mm*c(j));
                        end
                    end
                end
                u=u+d(j)*uu/sqrt(2*pi);
            end
        end
    case 'i' % inverse transform [data]=('i',sp,nincl,nazim)
        %or [data]=('a',sp,e,a) where e is a list of elevations and a is either a corresponding list of azimuths
        % or else a cell array contining several azimuths for each inclination
        % input data is sp=[(0,0) (1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) ... ]
        % length is (n+1)^2 where n is the order
        % output data is an array [ne,na]
        nsp=ceil(sqrt(length(a)))-1;
        [ue,we,lgp]=sphrharp(mv(3),b,nsp);
        if mv(3)=='a'
            na=nsp*2+1;
            ne=length(b);
        else
            na=c;
            ne=b;
        end
        i=0;
        da=zeros(ne,na);
        for nn=0:nsp
            % this could be vectorized somewhat to speed it up
            for mm=-nn:nn
                i=i+1;
                if i>length(a)
                    break;
                end
                da(:,1+mod(mm,na))=da(:,1+mod(mm,na))+a(i)*lgp(1+nn*(nn+1)/2+abs(mm),:)';
            end
        end
        if na>1 && mv(2)=='r' % convert real to complex only actually need to do this for min(b,floor(na/2)) values (but nyquist is a bit tricky)
            ix=2:ceil(na/2);
            iy=na:-1:na+2-ix(end);
            da(:,iy)=(da(:,ix)+1i*da(:,iy))*0.5;
            da(:,ix)=da(:,ix)-da(:,iy);
            da(:,1)=da(:,1)/sqrt(2);
        else
            da=da/sqrt(2);
        end
        if mv(3)=='a' % do a slow fft
            if iscell(c)
                u{ne,1}=[]; % reserve space for output cell array
                for i=1:ne
                    ai=c{i};
                    ui=repmat(da(i,1),1,length(ai));
                    for j=1:nsp
                        exj=exp(1i*j*ai(:).'); % we could vectorize this more
                        ui=ui+da(i,j+1).'.*exj+da(i,na+1-j).'.*conj(exj);
                    end
                    u{i}=ui/sqrt(pi);
                end
            else
                u=da(:,1);
                for j=1:nsp
                    exj=exp(1i*j*c(:)); % we could vectorize this more
                    u=u+da(:,j+1).*exj+da(:,na+1-j).*conj(exj);
                end
                u=u/sqrt(pi);
                if mv(2)=='r' && isreal(a)
                    u=real(u);
                end
            end
        else
            u=ifft(da,[],2)*na/sqrt(pi); % could put the scale factor 1/sqrt(pi) earlier
            if mv(2)=='r' && isreal(a)
                u=real(u);
            end
        end
    case 'c' % just output coordinates [inclination,azim,weights]=('c',nincl,nazim)
        % if m!='g', then order, a, must be even
        [u,w]=sphrharp(mv(3),a,0);
        if nargin<3
            b=ceil(a/1+(mv(3)=='g'));
        end
        v=(0:b-1)*2*pi/b;
end
if mv(4)~=' '
    switch mv(1)
        case 'f'
            if mv(2)=='r'
                ua=u;
                tit='Real Coefficients';
            else
                ua=abs(u);
                tit='Complex Coefficient Magnitudes';
            end
            nu=length(ua);
            no=ceil(sqrt(nu))-1;
            yi=zeros(no,2*no+1);
            for i=0:no
                for j=-i:i
                    iy=i^2+i+j+1;
                    if iy<=nu
                        yi(i+1,j+no+1)=ua(iy);
                    end
                end
            end
            imagesc(-no:no,0:no,yi);
            axis 'xy';
            if mv(4)=='P'
                colorbar;
            end
            xlabel('Harmonic');
            ylabel('Order');
            title(tit);
        case 'i'            % [data]=('i',sp,nincl,nazim)
            vv=(0:na)*2*pi/na;  % azimuth array
            gr=sin(ue)';
            surf(gr*cos(vv),gr*sin(vv),repmat(cos(ue)',1,na+1),real(u(:,[1:na 1])));
            axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            if mv(4)=='P'
                colorbar;
            end
            %             cblabel('Legendre Weight');
            %             title(sprintf('Sampling Grid: %s, %d, %d',mv(3),length(u),na-1));
        case 'c'            % [inclination,azim,weights]=('c',nincl,nazim)
            vv=[v v(1)];    % replicate initial azimuth point
            na=length(vv);
            gr=sin(u)';
            mesh(gr*cos(vv),gr*sin(vv),repmat(cos(u)',1,na),repmat(w',1,na));
            axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            if mv(4)=='P'
                colorbar;
                cblabel('Quadrature Weight');
            end
            title(sprintf('Sampling Grid: %s, %d, %d',mv(3),length(u),na-1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ueo,weo,lgpo]=sphrharp(gr,ne,nsp)
% calculate persistent variables for grid points and legendre polynomials
% we recalculate if transform order or inclination grid size or type are changed
% gr = grid type:
%       'u'=uniform inclination starting at pi/2n (default)
%       'U'=uniform but starting at north pole
%       'g'=gaussian inclination
%       'a'=arbitrary (specified by user - inverse transform only)
% ne = number of inclination values
% nsp = maximum order
% Outputs:
%    ueo(ne)     vector containing the inclination values
%    weo(ne)     vector containing the quadrature weights
%    lgpo((nsp+1)*(nsp+2)/2,ne)  evaluates the Schmidt seminormalized
%                                associated Legendre function at each value
%                                of cos(ue) for n=0:nsp, m=0:n.
persistent gr0 lgp ne0 ue we nsp0
if isempty(ne0)
    ne0=-1;
    nsp0=-1;
    gr0='a';
end
if gr=='a'
    ue=ne;
    ne=length(ue);
    ne0=-1;
    gr0=gr;
    nsp0=-1; % delete previous legendre polynomials
    lgp=[];
    we=[];
elseif gr~=gr0 || ne~=ne0
    if gr=='g'
        r = 1:ne-1;
        r = r ./ sqrt(4*r.^2 - 1);
        p = zeros( ne, ne );
        p( 2 : ne+1 : ne*(ne-1) ) = r;
        p( ne+1 : ne+1 : ne^2-1 ) = r;
        [q, ue] = eig(p);
        [ue, k] = sort(diag(ue));
        ue=acos(-ue)';
        we = 2*q(1,k).^2;
    elseif gr=='U'
        if rem(ne,2)
            error('inclination grid size must be even when using ''U'' option');
        end
        ue=(0:ne-1)*pi/ne;
        xx=zeros(1,ne);
        ah=ne/2;
        xx(1:ah)=(1:2:ne).^(-1);
        we=-4*sin(ue).*imag(fft(xx).*exp(-1i*ue))/ne;
    else % default is m='u'
        ue=(1:2:2*ne)*pi*0.5/ne;
        vq=(ne-abs(ne+1-2*(1:ne))).^(-1).*exp(-1i*(ue+0.5*pi));
        we=(-2*sin(ue).*real(fft(vq).*exp(-1i*(0:ne-1)*pi/ne))/ne);
    end
    gr0=gr;
    ne0=ne;
    nsp0=-1; % delete previous legendre polynomials
    lgp=[];
end
if nsp>nsp0
    lgp((nsp+1)*(nsp+2)/2,ne)=0; % reserve space
    for i=nsp0+1:nsp
        lgp(1+i*(i+1)/2:(i+1)*(i+2)/2,:)=legendre(i,cos(ue),'sch')*sqrt(0.5*i+0.25);
        lgp(1+i*(i+1)/2,:)=lgp(1+i*(i+1)/2,:)*sqrt(2);
    end
    nsp0=nsp;
end
lgpo=lgp;
weo=we;
ueo=ue;