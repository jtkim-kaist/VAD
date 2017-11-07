function [rgb,y,l]=v_colormap(map,m,n,p)
%V_COLORMAP set and plot color map
%
%   Usage: (1) v_colormap([],'g'); % plot the current color map
%          (2) v_colormap([],'',256); % intepolate the current color map to 256 entries
%          (3) v_colormap('copper','y'); % make copper color map linear in luminance
%          (4) v_colormap('copper','L'); % make copper color map linear in lightness^2
%
%  Inputs:
%           map  Either an (r,3) array specifying the RGB colourmap entries
%                or else a string specifying a bult-in colour map. Use []
%                to leave the colour map unchanged.
%                Standard maps are:
%                   'Jet','HSV','Hot','Cool','Spring','Summer','Autumn',
%                   'Winter','Gray','Bone','Copper','Pink','Lines','Parula'
%                Additional maps, all with 64 entries, are:
%                   'v_thermliny'  thermal scale that is linear in luminance
%                                  Varies from black through blue, red, yellow to white.
%                                  Luminance varies linearly from black to white.
%                   'v_bipliny'    bipolar scale that is linear in luminance
%                                  Negative values are blue/turquoise and postive value are orange/yellow.
%                                  Luminance varies linearly from black to white with zero at 50.8% gray.
%                   'v_bipveey'    bipolar scale that is V-shaped in luminance
%                                  Negative values are blue/turqoise and positive values are red/yellow.
%                                  Luminance is proportional to absolute value with zero=black.
%                For the two bipolar scales, zero corresponds to entry 33 and so the range of values
%                is -32:31 or, equivalently, either -1 to +0.96875 or -1.0323 to +1.
%
%             m  Mode string:
%                   'g' to plot information about the color map
%                   'y' to force luminance^p to be linear or V-shaped (two linear segments)
%                   'l' to force lightness^p to be linear or V-shaped (two linear segments)
%                   'Y' like 'y' but with default p=0.667
%                   'L' like 'l' but with default p=2
%                   'f' flips the map to reverse its order
%                   'b' make maximum luminance 0.05 (or 0.1 for 'B')
%                   'w' make maximum luminance 0.95 (or 0.9 for 'W')
%
%             n  the number of entries in the colourmap or the number in
%                each linearly-interpolated segment excluding the entry shared
%                with the previous segment. The total number of entries is n=sum(n).
%                For modes 'y','Y','l','L' the number of segments must be 1
%                or 2; otherwise the number of segments must be 1 or r-1.
%
%             p  power law to use for linearized luminance or lightness [default p=1]
%                see the description of 'y' and 'l' for its effect
%
% Outputs:
%           rgb  RGB color map entries; one per row.
%                All values will be in the range 0 to 1
%
%             y  column vector of luminance values
%
%             l  column vector of lightness values (lightness is the perceived brightness)

% Bugs/Suggestions:
% (1) add option to exclude black from the colormap

%      Copyright (C) Mike Brookes 2012
%      Version: $Id: v_colormap.m 7758 2016-04-12 07:11:33Z dmb $
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
persistent maps nams mcal modes nszs pows la lb lc lci laci lk lq yv
if isempty(maps)
    lk=(6/29)^3;
    la=841/108;
    lb=4/29;
    lc=1.16;
    lq=la*lc*lk;
    lci=1/lc;
    laci=lci/la;
    %     yv=inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);
    %     yv=yv(:,1);
    yv=[0.2126; 0.7152; 0.0722];
    nams={'v_thermliny','v_bipliny','v_bipveey'};
    % modified thermal with better grayscale linearity
    mcal=[1 1 1]; % +1 if need to calculate maps entry
    modes={'y' 'y' 'y'}; % modes for map calculation
    nszs={64 64 [33 31]};  % sizes for maps
    pows=[1 1 1];  % powers for maps
    maps{1}=[0 0 0; 72 0 167; 252 83 16; 255 249 0; 255 255 255]/255;
    % bipolar map with grayscale linearity
    maps{2}=[0 0 0; 0 2 203; 1 101 226; 128 128 128; 252 153 12; 252 245 0; 252 249 18; 252 252 252]/252;
    maps{3}=[0 0.95 1; 0 0 0.9; 0 0 0; 0.5 0 0; 0.80 0.78 0];
end
if nargin<2
    m='';
end
if nargin<4
    if any(m=='Y')
        p=2/3;
    elseif any(m=='L')
        p=2;
    else
        p=1;  % power to raise lightness/luminance to
    end
end
pr=1/p;
um=m;
m=lower(m);   % convert mode letters to lower case
gotmap= nargin==0 || numel(map)==0;   % already got the map
if ~gotmap && ischar(map)   % if map given as a string
    ix=find(strcmpi(map,nams),1); % check if it is one of ours
    if numel(ix)
        if mcal(ix)  % need to calculate the map the first time around
            maps{ix}=v_colormap(maps{ix},modes{ix},nszs{ix},pows(ix));
            mcal(ix)=0;  % don't calculate it again
        end
        colormap(maps{ix});
    else
        colormap(map); % not one of ours - just pass it on to standard colormap function
    end
    gotmap=1;
end
if any(m=='y') ||  any(m=='l') || (nargin>2 && numel(n)>0) % need to do linear interpolation
    if gotmap
        map=colormap;  % linearly interpolate the existing map
    end
    nm=size(map,1);
    if any(m=='y') ||  any(m=='l')
        y=map*yv;  % convert map to luminance
        up=y(2:nm)>y(1:nm-1); % find increasing
        ex=up(1:nm-2)~=up(2:nm-1); % +1 for a peak or valley
        yd=2*up(1)-1; % +1 if y initially increasing
        switch sum(ex)
            case 0 % monotonic
                if nargin<3
                    r=nm;
                else
                    r=n(1);
                end
                if any(m=='y')
                    l=y([1 nm]).^p;
                    tt=(l(1)+(0:r-1)'*(l(2)-l(1))/(r-1)).^pr; % target luminances
                else
                    tt=y([1 nm]');
                    l=(lc*(la*tt+(tt>lk).*(tt.^(1/3)-la*tt-lb))).^p;
                    tt=(l(1)+(0:r-1)'*(l(2)-l(1))/(r-1)).^pr; % target lightnesses
                    tt=laci*tt+(tt>lq).*((lci*tt+lb).^3-laci*tt); % target luminances
                end
                [ss,ix]=sort([tt;y]*yd);
            case 1 % V-shaped
                ipk=find(ex,1)+1;  % position of peak/valley in y
                if nargin<3
                    n=[ipk nm-ipk];  % size of linear segments
                end
                r=n(1)+n(2);  % total size of map
                if any(m=='y')
                    l=y([1 ipk nm]).^p;
                    tt=(l(2)+[(1-n(1):0)*(l(2)-l(1))/(n(1)-1) (1:n(2))*(l(3)-l(2))/(n(2))]').^pr; % target luminances
                else
                    tt=y([1 ipk nm]');
                    l=(lc*(la*tt+(tt>lk).*(tt.^(1/3)-la*tt-lb))).^p;
                    tt=(l(2)+[(1-n(1):0)*(l(2)-l(1))/(n(1)-1) (1:n(2))*(l(3)-l(2))/(n(2))]').^pr; % target lightnesses
                    tt=laci*tt+(tt>lq).*((lci*tt+lb).^3-laci*tt); % target luminances
                end
                [ss,ix]=sort([tt(1:n(1))-y(ipk); y(ipk)-tt(n(1)+1:r);y(1:ipk)-y(ipk); y(ipk)-y(ipk+1:nm)]*yd);
            otherwise
                error('luminance has more than two monotonic segments');
        end
    else                % just linearly interpolate the given values
        if numel(n)==nm-1
            r=sum(n);
            y=[1;cumsum(n(:))];
        else
            r=n(1);
            y=1+(0:nm-1)'*(r-1)/(nm-1);
        end
        tt=(1:r)';
        [ss,ix]=sort([tt;y]);
    end
    jx=zeros(size(ix));
    jx(ix)=1:numel(jx);
    jx=min(max(jx(1:r)-(1:r)',1),nm-1);
    al=(tt-y(jx))./(y(jx+1)-y(jx)); % fraction of upper sample to include
    map=map(jx,:)+(map(jx+1,:)-map(jx,:)).*al(:,ones(1,3));
    colormap(map);
    gotmap=1;
end
if ~gotmap
    colormap(map);  % directly specified numerical map
end
if any(m=='f')
    rgb=colormap; % get output values from the current colourmap
    colormap(rgb(end:-1:1,:))
end
rgb=colormap; % get output values from the current colourmap
y=rgb*yv;  % convert RGB to luminance
minyt=0.05*(any(m=='b')+any(um=='B')); % target minimum luminance
maxyt=1-0.05*(any(m=='w')+any(um=='W')); % target maximum luminance
maxy=max(y);
miny=min(y);
if maxy>maxyt || miny<minyt
    maxy=max(maxy,maxyt);
    miny=min(miny,minyt);
    rgb=(rgb-miny)*(maxyt-minyt)/(maxy-miny)+minyt;
    colormap(rgb);
    y=rgb*yv;  % convert RGB to luminance
end
l=lc*(la*y+(y>lk).*(y.^(1/3)-la*y-lb)); % convert luminance to lightness
if any(m=='g')
    sp=[1 2 2];
    ssp=sum(sp);
    axw=0.05;
    nc=size(rgb,1);  % size of color map
    hsv=rgb2hsv(rgb);
    subplot(ssp,1,sp(1)+(1:sp(2)));
    plot(1:nc,y,'-k');
    hold on
    plot(1:nc,rgb(:,1),'-r');
    plot(1:nc,rgb(:,2),'-g');
    plot(1:nc,rgb(:,3),'-b');
    hold off
    axis([0.5 nc+0.5 -axw 1+axw]);
    ylabel('RGB + Y');
    subplot(ssp,1,sp(1)+sp(2)+1:ssp);
    plot(1:nc,l,'-k');
    hold on
    plot(1:nc,hsv(:,1),'-r');
    plot(1:nc,hsv(:,2),'-g');
    plot(1:nc,hsv(:,3),'-b');
    hold off
    axis([0.5 nc+0.5 -axw 1+axw]);
    ylabel('HSV + L*');
    subplot(ssp,1,1:sp(1));
    image(permute(reshape([rgb y(:,[1 1 1])],[nc,3,2]),[3 1 2]));
    set(gca,'YTick',[]);
end
