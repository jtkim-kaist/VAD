function [gci goi] = v_sigma(lx,fs,fmax)

%   Singularity in EGG by Multiscale Analysis (SIGMA) Algorithm
%
%   [gci goi] = v_sigma(lx,fs,fmax)
%
%   Inputs:
%       lx      Nx1 vector LX signal
%       fs      Sampling freq (Hz)
%       fmax    [Optional] max laryngeal freq
%   Outputs:
%       gci      Vector of gcis as sample instants.
%       goi      Vector of gois as sample instants.
%
%   External Functions:
%       Function gaussmix in Voicebox.
%
%   Notes:
%       Due to the initialization of the EM algorithm on a random data
%       point, V_SIGMA does not always produce deterministic behaviour.
%
%   References:
%       M. R. P. Thomas and P. A. Naylor, "The SIGMA Algorithm: A Glottal
%       Activity Detector for Electroglottographic Signals," IEEE Trans.
%       Audio, Speech, Lang. Process., vol.17, no.8, pp.1557-1566, Nov.
%       2009.
%
%   Revision History:
%       13/09/10: Swallow postfilter threw an error when fewer than 3 GCIs 
%       were detected. Final GOI cycle threw an error when search bounds
%       exceeded length of signal. Similar problem in line 179 fixed.
%
%**************************************************************************
% Author:           M. R. P. Thomas, 29 May 2007
%**************************************************************************
%      Copyright (C) M. R. P. Thomas, Mike Brookes 2007-2013
%      Version: $Id: v_sigma.m 4752 2014-06-25 12:56:14Z dmb $
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

if(nargin<3)
    fmax = 400;
end

% PARAMETERS - we don't like these
fmin = 50;      % GOI post-filtering
Tmin = 1/fmax;  % Sets GD window length
Tmax = 1/fmin;  % GOI post-filtering
oqmin = 0.1;    % GOI post-filtering
oqmax = 0.9;    % GOI post-filtering

gwlen = Tmin;   % group delay evaluation window length (normally 0.0025)
fwlen=0.000;    % window length used to smooth group delay (normally 0.000)

% Normalise (for plotting purposes)
lx = lx/max(abs(lx));

% Calculate SWT
nlev = 3;
nu = length(lx);
nU=(2^nlev)*ceil(nu./(2^nlev));
[Lo_D Hi_D] = wfilters('bior1.5','d');
[swa swd] = swt([lx; zeros(nU-nu,1)],nlev,Lo_D,Hi_D);
swa=swa(:,1:nu); swd=swd(:,1:nu);
mp = prod(swd)';
mp = [0;mp];

% Find third roots
nmp = mp;
nmp(find(nmp>0)) = 0;   % Half-wave rectify on negative half of mp for GCI
crnmp = nthroot(nmp,3);

pmp = mp;
pmp(find(pmp<0)) = 0;   % Half-wave rectify on positive half of mp for GOI
crpmp = nthroot(pmp,3);

% Group delay evaluation on mp
[gcic sew ngrdel ntoff] = xewgrdel(nmp,fs,gwlen,fwlen);
ngrdel=-[zeros(ntoff,1); ngrdel(1:end-ntoff)];
[goic sew pgrdel ptoff] = xewgrdel(pmp,fs,gwlen,fwlen);
pgrdel=-[zeros(ptoff,1); pgrdel(1:end-ptoff)];

% Set up other variables
gci = zeros(1,length(lx));      
gciall = zeros(1,length(lx));
gciall(:) = NaN;
gciall(round(gcic)) = 0.25;

goi = zeros(1,length(lx));      
goiall = zeros(1,length(lx));
goiall(:) = NaN;
goiall(round(goic)) = 0.25;

% --------- GCI Detection ---------- %

% Model GD slope
mngrdel = (-((gwlen*fs)/2-1):(gwlen*fs)/2-1)';
mngrdellen = length(mngrdel);
cmngrdel = zeros(length(ngrdel),1);

for i=1:length(gcic)
    lbnd = round(gcic(i)-((gwlen*fs)/2-1)); 
    ubnd = lbnd+mngrdellen-1;
    
    if ~( (lbnd<1) || (ubnd>length(ngrdel)) )
        nfv(i,1) = sum(crnmp(lbnd:ubnd));                     % Sum of crnmp over GD window
        nfv(i,2) = min(crnmp(lbnd:ubnd));                     % Peak value of crnmp       
        nfv(i,3) = sqrt( mean( (mngrdel-ngrdel(lbnd:ubnd)).^2 ) ); % Phase slope deviation        
        
        cmngrdel(lbnd:ubnd) = mngrdel;
    end
end

nclust = 2;
snfv = size(nfv);

% Determine clusters
[mm,vv,ww]=gaussmix(nfv,var(nfv)/100,[],nclust,'kf');

% Find cluster with lowest crnmp sum
[y I] = min(mm(:,1));

% Find log likelihoods for each cluster
for i=1:nclust
    [m_,v_,w_,g,f,ll(i,:)]=gaussmix(nfv,[],0,mm(i,:),vv(i,:),ww(i,:));
end
[m,in]=max(ll); % Find which cluster each feature vector belongs to

% close all; % commented out by dmb, 12/11/2013
naccept = [];
nreject = [];
for i=1:snfv(1)
    if (in(i) == I)
        naccept = [naccept; nfv(i,:)];
    else
        nreject = [nreject; nfv(i,:)];
    end
end

% If the candidate belongs to the chosen cluster then keep
for i=1:length(nfv)
    if (in(i) == I)
        gci(round(gcic(i))) = 0.5;
    end
end

% -------- Post-filter swallows (GCIs only) -------- %
if(length(gci)>2)
    % If a gci is separated from all others by more than Tmax, delete.
    fgci = find(gci);
    % Check first one
    if ( (fgci(2)-fgci(1)) > Tmax*fs)
        fgci = fgci(2:end);
    end
    % Check the middle
    i=2;
    while(i<length(fgci)-1)
        if ( ((fgci(i)-fgci(i-1))>Tmax*fs) && ((fgci(i+1)-fgci(i))>Tmax*fs) )
            fgci = [fgci(1:i-1) fgci(i+1:end)];
        end
        i = i+1;
    end
    % Check last one
    if ( (fgci(end)-fgci(end-1)) > Tmax*fs)
        fgci = fgci(1:end-1);
    end
    % Convert back
    gci = zeros(1,max(fgci));
    gci(fgci) = 0.5;
end

% --------- GOI Detection ---------- %
% Model GD slope
mpgrdel = (-((gwlen*fs)/2-1):(gwlen*fs)/2-1)';
mpgrdellen = length(mpgrdel);
cmpgrdel = zeros(length(pgrdel),1);

for i=1:length(goic)
    lbnd = round(goic(i)-((gwlen*fs)/2-1)); 
    %ubnd = round(goic(i)+((gwlen*fs)/2-1)); 
    ubnd = min(lbnd+mpgrdellen-1,length(crpmp));
    
    if ~( (lbnd<1) || (ubnd>length(pgrdel)) )
        pfv(i,1) = sum(crpmp(lbnd:ubnd));                     % Sum of crnmp over GD window
        pfv(i,2) = max(crpmp(lbnd:ubnd));                     % Peak value of crpmp
        pfv(i,3) = sqrt( mean( (mpgrdel-pgrdel(lbnd:ubnd)).^2 ) ); % Phase slope deviation
                
        cmpgrdel(lbnd:ubnd) = mpgrdel;
    end
end

nclust = 2;
spfv = size(pfv);

% Determine clusters
[mm,vv,ww]=gaussmix(pfv,var(pfv)/100,[],nclust,'kf');

% Find cluster with highest crpmp sum
[y I] = max(mm(:,1));

% Find log likelihoods for each cluster
ll = [];
for i=1:nclust
    [m_,v_,w_,g,f,ll(i,:)]=gaussmix(pfv,[],0,mm(i,:),vv(i,:),ww(i,:));
end
[m,in]=max(ll); % Find which cluster each feature vector belongs to

paccept = [];
preject = [];
for i=1:spfv(1)
    if (in(i) == I)
        paccept = [paccept; pfv(i,:)];
    else
        preject = [preject; pfv(i,:)];
    end
end

% If the candidate belongs to the chosen cluster then keep
for i=1:length(pfv)
    if (in(i) == I)
        goi(round(goic(i))) = 0.75;
    end
end

% ------- GOI Post-Filtering ------- %

% For all GCIs, find GOIs which are within OQ limits
goiprefilt = goi;
goifilt = zeros(size(goi));
gciind = find(gci);
Tprev = Tmax;
prev = 0;
nofirst = 0;
for i=2:length(gciind)
    lbnd = gciind(i-1);
    ubnd = gciind(i);
    T = ubnd-lbnd;
    if(T>Tmax*fs)   % If period is too long, use previous.
        T = Tprev;
    end
    I = find(goi( round(lbnd+(1-oqmax)*T):round(lbnd+(1-oqmin)*T) ));
    if(~isempty(I)) % If we have a GOI
        prev = round(I(1)+(1-oqmax)*T-1);
        goifilt(round(I(1)+lbnd+(1-oqmax)*T-1)) = 0.5;  % Taking first - should it be highest energy?
    else            % If not then insert at last OQ.
        if(prev>0)
            if( (lbnd+prev) < gciind(i) ) % Protect against GOI occuring after next GCI
                goifilt(lbnd + prev-1) = 0.5;
            else
                goifilt( gciind(i)-1) = 0.5;
            end
        end
        if(i==2)
            nofirst = 1;
        end
    end
    if(nofirst && (prev>0)) % If no GOI occurs after GOI, no previous period exists, so add after a period has been found.
        goifilt(gciind(1)+prev-1) = 0.5;
        nofirst = 0;
    end
    Tprev = T;
end
% Final period
lbnd = gciind(end);
I = find(goi( round(lbnd+(1-oqmax)*T):min(round(lbnd+(1-oqmin)*T),length(goi))));
if(~isempty(I))
    goifilt(round(I(1)+lbnd+(1-oqmax)*T-1)) = 0.5;
else            % If not then insert at last OQ.
    if(prev>0)
        goifilt(lbnd + prev) = 0.5;
    end
end
goi = goifilt;

gci = find(gci>0);
goi = find(goi>0);

%% EW group delay epoch extraction
function [tew,sew,y,toff]=xewgrdel(u,fs,gwlen,fwlen)

error(nargchk(2,4,nargin));

if(nargin < 4)
    fwlen = voicebox('dy_fwlen');          % window length to smooth group delay
end
if (nargin < 3)
    gwlen = voicebox('dy_gwlen');          % window length of group delay
end

% perform group delay calculation

gw=2*floor(gwlen*fs/2)+1;            % force window length to be odd
ghw=window('hamming',gw,'s');
ghw = ghw(:);                           % force to be a column (dmb thinks window gives a row - and he should know as he wrote it!)
ghwn=ghw'.*(gw-1:-2:1-gw)/2;            % weighted window: zero in middle

u2=u.^2;
yn=filter(ghwn,1,u2);
yd=filter(ghw,1,u2);
yd(abs(yd)<eps)=10*eps;                 % prevent infinities
y=yn(gw:end)./yd(gw:end);               % delete filter startup transient
toff=(gw-1)/2;
fw=2*floor(fwlen*fs/2)+1;            % force window length to be odd
if fw>1
    daw=window('hamming',fw,'s');
    y=filter(daw,1,y)/sum(daw);         % low pass filter 
    toff=toff-(fw-1)/2;
end
[tew,sew]=zerocros(y,'n');              % find zero crossings

tew=tew+toff;                           % compensate for filter del