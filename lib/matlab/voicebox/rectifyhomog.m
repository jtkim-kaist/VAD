function [imr,xa,ya]=rectifyhomog(ims,roc,k0,mode)
%RECTIFYHOMOG Apply rectifying homographies to an image set
%
% Usage:    figure(101);     % Initial figure for rectiied image display
%           rectifyhomog(ims,roc,k0,'ga');   % plot in individual figures
%
% Inputs:
%        ims{nc}       cell array of input images (colour or monochrome)
%        roc(3,3,nc)   rotation matrices from world coordinates to camera coordinates
%        k0            camera matrix or focal length in pixels optionally divided by the image width [0.8]
%        mode          mode string
%                         g  show images on separate figures
%                         G  tile images onto a single figure [default if no output arguments]
%                         k  clip to original image dimensions
%                         l  do not link axes
%                         v  k0 contains diffent values for each camera
%                         a  orient to average camera orientation
%                        [z  make optical axis perpendicular to camera plane]
%                        [x  align x axis with camera displacements]
% Outputs:
%        imr{nc}(my,mx,nc)  output images (uint8)
%        xa{nc}(mx)        x axis for each image
%        ya{nc}(my)        y axis for each image

%      Copyright (C) Mike Brookes 2012
%      Version: $Id: rectifyhomog.m 1642 2012-03-16 16:22:26Z dmb $
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
if ~iscell(ims)
    ims={ims};
end
nc=numel(ims); % number of images
if ~((nc==1 && isequal(size(roc),[3 3])) || isequal(size(roc),[3 3 nc]))
    error('roc should have dimensions [3,3,%d]',nc);
end
if nargin<4 || ~numel(mode)
    mode='';
end
if nargin<3 || ~numel(k0)
    k0=0.8;
end
if nargin<2 || ~numel(roc)
    roc=repmat(eye(3),[1 1 nc]);
end

% sort out options

vv=any(mode=='v');  % different k0 for each camera
nk=numel(k0)/(1+vv*(nc-1));  % number of k parameters per camera
switch nk
    case {1,3}
        k0=reshape(k0,[1 nk 1+vv*(nc-1)]);
    case 9
        k0=reshape(k0,[3 3 1+vv*(nc-1)]);
    otherwise
        error('k0 must have 1, 3 or 9 values per camera');
end
ncr=1+(nargout>0)*(nc-1);   % numer of output images
if any(mode=='g')
    gmode=1;
elseif any(mode=='G') || ~nargout
    gmode=2;
else
    gmode=0;
end
if gmode>0
    fig0=gcf; % initialize the figure
end
if any(mode=='k')
    modeh='kxt';
else
    modeh='t';
end

% determine a global camera rotation

if any(mode=='a')
    qrc=zeros(4,nc);  % calculate the mean camera orientation
    for i=1:nc
        qrc(:,i)=rotro2qr(roc(:,:,i));
    end
    rocmean=rotqr2ro(rotqrmean(qrc));
else
    rocmean=eye(3);
end

% now do image transformations

imr=cell(ncr,1);
xa=imr;
ya=imr;
axh=zeros(nc,1);
splx=ceil(sqrt(nc));
sply=ceil(nc/splx);
for i=1:nc
    k0i=k0(:,:,1+vv*(i-1));  % camera parameters
    if numel(k0i)<9
        imsz=size(ims{i});
        fe=k0(1);
        if fe<0.1*imsz(2) % focal length is a fraction of the width
            fe=k0*imsz(2);
        end
        if numel(k0i)<3
            xy0=(imsz(2:-1:1)+1)/2;
        else
            xy0=k0i(2:3);
        end
        k0=eye(3);
        k0(1:4:5)=fe;
        k0(7:8)=xy0;
    end
    j=min(i,ncr);
    rocall=rocmean*roc(:,:,i)';
    titl=sprintf('%d: pan-tilt-roll = %.1f°, %.1f°, %.1f°',i,-rotro2eu('yxz',rocall)*180/pi);
    [imr{j},xa{j},ya{j}]=imagehomog(uint8(ims{i}),k0*rocall/k0,modeh);  % apply inverse of rotation matrix
    if gmode>0
        if gmode>1
            subplot(sply,splx,i);
        else
            figure(fig0+i-1);
        end
        imagesc(xa{j},ya{j},imr{j});
        axis image
        title(titl);
        axh(i)=gca;
    end
end
if gmode>0 && nc>1
    if ~any(mode=='l')
        linkaxes(axh);
    end
    figure(fig0);
end
