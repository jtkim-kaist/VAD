function [gci,goi] = dypsa(s,fs)
%DYPSA   Derive glottal closure instances from speech [gci,goi] = (s,fs)
%   Note: Needs to be combined with a voiced-voiceless detector to eliminate
%   spurious closures in unvoiced and silent regions.
%
%   Inputs:
%   s     is the speech signal
%   fs    is the sampling frequncy
%
%   Outputs:
%   gci   is a vector of glottal closure sample numbers
%   gco   is a vector of glottal opening sample numbers derived from
%         an assumed constant closed-phase fraction
%
%   References:
%   [1]  P. A. Naylor, A. Kounoudes, J. Gudnason, and M. Brookes, “Estimation of Glottal Closure
%        Instants in Voiced Speech using the DYPSA Algorithm,” IEEE Trans on Speech and Audio
%        Processing, vol. 15, pp. 34–43, Jan. 2007.
%   [2]  M. Brookes, P. A. Naylor, and J. Gudnason, “A Quantitative Assessment of Group Delay Methods
%        for Identifying Glottal Closures in Voiced Speech,” IEEE Trans on Speech & Audio Processing,
%        vol. 14, no. 2, pp. 456–466, Mar. 2006.
%   [3]  A. Kounoudes, P. A. Naylor, and M. Brookes, “The DYPSA algorithm for estimation of glottal
%        closure instants in voiced speech,” in Proc ICASSP 2002, vol. 1, Orlando, 2002, pp. 349–352.
%   [4]  C. Ma, Y. Kamp, and L. F. Willems, “A Frobenius norm approach to glottal closure detection
%        from the speech signal,” IEEE Trans. Speech Audio Processing, vol. 2, pp. 258–265, Apr. 1994.
%   [5]  A. Kounoudes, “Epoch Estimation for Closed-Phase Analysis of Speech,” PhD Thesis,
%        Imperial College, 2001.

% Algorithm Parameters
%       The following parameters are defined in voicebox()
%
%   dy_cpfrac=0.3;           % presumed closed phase fraction of larynx cycle
%   dy_cproj=0.2;            % cost of projected candidate
%   dy_cspurt=-0.45;         % cost of a talkspurt
%   dy_dopsp=1;              % Use phase slope projection (1) or not (0)?
%   dy_ewdly=0.0008;         % window delay for energy cost function term [~ energy peak delay from closure] (sec)
%   dy_ewlen=0.003;          % window length for energy cost function term (sec)
%   dy_ewtaper=0.001;        % taper length for energy cost function window (sec)
%   dy_fwlen=0.00045;        % window length used to smooth group delay (sec)
%   dy_fxmax=500;            % max larynx frequency (Hz) 
%   dy_fxmin=50;             % min larynx frequency (Hz) 
%   dy_fxminf=60;            % min larynx frequency (Hz) [used for Frobenius norm only]
%   dy_gwlen=0.0030;         % group delay evaluation window length (sec)
%   dy_lpcdur=0.020;         % lpc analysis frame length (sec)
%   dy_lpcn=2;               % lpc additional poles
%   dy_lpcnf=0.001;          % lpc poles per Hz (1/Hz)
%   dy_lpcstep=0.010;        % lpc analysis step (sec)
%   dy_nbest=5;              % Number of NBest paths to keep
%   dy_preemph=50;           % pre-emphasis filter frequency (Hz) (to avoid preemphasis, make this very large)
%   dy_spitch=0.2;           % scale factor for pitch deviation cost
%   dy_wener=0.3;            % DP energy weighting
%   dy_wpitch=0.5;           % DP pitch weighting
%   dy_wslope=0.1;           % DP group delay slope weighting
%   dy_wxcorr=0.8;           % DP cross correlation weighting
%   dy_xwlen=0.01;           % cross-correlation length for waveform similarity (sec)

%   Revision History: 
%
%   3.0 - 29 Jun 2006  - Rewrote DP function to improve speed
%   2.6 - 29 Jun 2006  - Tidied up algorithm parameters
%   2.4 - 10 Jun 2006  - Made into a single file aand put into VOICEBOX
%   2.3 - 18 Mar 2005  - Removed 4kHz filtering of phase-slope function 
%   2.2 - 05 Oct 2004  -  dpgci uses the slopes returned from xewgrdel
%                      -  gdwav from speech with fs<9000 is not filtered
%                      -  Various outputs and inputs of functions have been
%                         removed since now there is no plotting
%   1.0 - 30 Jan 2001  - Initial version [5]

%   Bugs:
%         1. Allow the projections only to extend to the end of the larynx cycle
%         2. Compensate for false pitch period cost at the start of a voicespurt
%         3. Should include energy and phase-slope costs for the first closure of a voicespurt
%         4. should delete candidates that are too close to the beginning or end of speech for the cost measures
%            currently this is 200 samples fixed in the main routine but it should adapt to window lengths of
%            cross-correlation, lpc and energy measures.
%         5. Should have an integrated voiced/voiceless detector
%         6. Allow dypsa to be called in chunks for a long speech file
%         7. Do forward & backward passes to allow for gradient descent and/or discriminative training
%         8. Glottal opening approximation does not really belong in this function
%         9. The cross correlation window is asymmetric (and overcomplex) for no very good reason
%        10. Incorporate -0.5 factor into dy_wxcorr and abolish erroneous (nx2-1)/(nx2-2) factor
%        11. Add in talkspurt cost at the beginning rather than the end of a spurt (more efficient)
%        12. Remove qmin>2 condition from voicespurt start detection (DYPSA 2 compatibility) in two places
%        13. Include energy and phase-slope costs at the start of a voicespurt
%        14. Single-closure voicespurt are only allowed if nbest=1 (should always be forbidden)
%        15. Penultimate closure candidate is always acceptd
%        16. Final element of gcic, Cfn and Ch is unused
%        17. Needs to cope better with irregular voicing (e.g. creaky voice)
%        18. Should give non-integer GCI positions for improved accuracy
%        19. Remove constraint that first voicespurt cannot begin until qrmax after the first candidate

%      Copyright (C) Tasos Kounoudes, Jon Gudnason, Patrick Naylor and Mike Brookes 2006
%      Version: $Id: dypsa.m 3102 2013-06-13 21:16:11Z dmb $
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

% Extract algorithm constants from VOICEBOX


dy_preemph=voicebox('dy_preemph');
dy_lpcstep=voicebox('dy_lpcstep');
dy_lpcdur=voicebox('dy_lpcdur');
dy_dopsp=voicebox('dy_dopsp');              % Use phase slope projection (1) or not (0)?
dy_ewtaper=voicebox('dy_ewtaper');        % Prediction order of FrobNorm method  in seconds
dy_ewlen=voicebox('dy_ewlen');        % windowlength of FrobNorm method  in seconds
dy_ewdly=voicebox('dy_ewdly');        % shift for assymetric speech shape at start of voiced cycle
dy_cpfrac=voicebox('dy_cpfrac');        % presumed ratio of larynx cycle that is closed
dy_lpcnf=voicebox('dy_lpcnf');          % lpc poles per Hz (1/Hz)
dy_lpcn=voicebox('dy_lpcn');            % lpc additional poles
dy_xwlen=voicebox('dy_xwlen');            % cross-correlation length for waveform similarity (sec)
dy_fxminf=voicebox('dy_fxminf');            % minimum pitch for Frobenius norm calculation

lpcord=ceil(fs*dy_lpcnf+dy_lpcn);       % lpc poles

%PreEmphasise input speech
s_used=filter([1 -exp(-2*pi*dy_preemph/fs)],1,s);

% perform LPC analysis, AC method with Hamming windowing
[ar, e, k] = lpcauto(s_used,lpcord,floor([dy_lpcstep dy_lpcdur]*fs));

if any(any(isinf(ar)))    % if the data is bad and gives infinite prediction coefficients we return with a warning
    warning('No GCIs returned');
    gci=[];
    return;
end;

% compute the prediction residual
r = lpcifilt(s_used,ar,k); 

% compute the group delay function:  EW method from reference [2] above
[zcr_cand,sew,gdwav,toff]=xewgrdel(r,fs); 
gdwav=-[zeros(toff,1); gdwav(1:end-toff)];
zcr_cand=[round(zcr_cand), ones(size(zcr_cand))];   %flag zero crossing candidates with ones

sew=0.5+sew';  %the phase slope cost of each candidate

pro_cand=[];
if dy_dopsp ~= 0
    pro_cand = psp(gdwav,fs);
    pro_cand = [pro_cand, zeros(length(pro_cand),1)]; %flag projected candidates with zeros
    sew =      [sew zeros(1,size(pro_cand,1))];      %the phase slope cost of a projected candidate is zero
end;

%Sort the zero crossing and projected candidates together and remove any candidates that
%are to close to the start or end of the speech signal because the cost functions
%need room either side. 

[gcic,sin] = sortrows([zcr_cand; pro_cand],1);  
sew=sew(sin);
saf=max([200,dy_xwlen*fs/2+1,fs/dy_fxminf]);
sin=find(and(saf<gcic,gcic<length(gdwav)-saf));
gcic=gcic(sin,:);
sew=sew(sin);

% compute the frobenious norm function used for a cost in the DP
fnwav=frobfun(s_used,dy_ewtaper*fs,dy_ewlen*fs,dy_ewdly*fs);

%Dynamic programming, picks the most likely candidates based on the pitch consistency, energy etc.
[gci] = dpgci(gcic, s_used(:), sew, fnwav, fs);

%Evaluate goi ... determined to be dy_cpfrac percentage of the larynx cylce away from the last gci
goi=simplegci2goi(gci,dy_cpfrac);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = psp(g,fs)
%PSP  Calculates the phase slope projections of the group delay function
%   Z = PSP(G) computes the 

%   Revision History: 
%       V1.0 July 12th 2002:
%            Nov. 6th 2002 : if statement added to remove "last" midpoint 

g = g(:);

gdot = [diff(g);0];
gdotdot = [diff(gdot);0];

% find the turning points  as follows: [tp_number, index_of_tp, min(1) or max(-1), g(index_of_tp)]
turningPoints = zcr(gdot);
turningPoints = [[1:length(turningPoints)]', turningPoints, sign(gdotdot(turningPoints)), g(turningPoints)];

% useful for debug/plotting
%tplot = zeros(length(g),1);
%tplot(turningPoints(:,1)) = turningPoints(:,2);

% find any maxima which are < 0 
negmaxima = turningPoints(find(turningPoints(:,3) == -1 & turningPoints(:,4) < 0 & turningPoints(:,1)~=1),:);  %Change 01.05.2003 JG: The first row can't be included

% find the midpoint between the preceding min and the negative max
nmi = negmaxima(:,1);
midPointIndex = turningPoints(nmi-1,2) + round(0.5*(turningPoints(nmi,2) - turningPoints(nmi-1,2)));
midPointValue = g(midPointIndex);

% project a zero crossing with unit slope
nz = midPointIndex - round(midPointValue);

% find any minima which are > 0 
posminima = turningPoints(find(turningPoints(:,3) == 1 & turningPoints(:,4) > 0),:);

% find the midpoint between the positive min and the following max
pmi = posminima(:,1); 

%Remove last midpoint if it is the last sample
if ~isempty(pmi), if pmi(end)==size(turningPoints,1), pmi=pmi(1:end-1); end; end;

midPointIndex = turningPoints(pmi,2) + round(0.5*(turningPoints(pmi+1,2) - turningPoints(pmi,2)));
midPointValue = g(midPointIndex);

% project a zero crossing with unit slope
pz = midPointIndex - round(midPointValue);

z = sort([nz;pz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function i = zcr(x, p)
%ZCR  Finds the indices in a vector to  zero crossings
%   I = ZCR(X) finds the indices of vector X which are closest to zero-crossings.
%   I = ZCR(X, P) finds indices for positive-going zeros-crossings for P=1 and
%   negative-going zero-crossings for P=0.

x = x(:);

if (nargin==2)
    if (p==0) 
        z1 = zcrp(x);   % find positive going zero-crossings
    elseif (p==1) 
        z1 = zcrp(-x);  % find negative going zero-crossings
    else
        error('ZCR: invalid input parameter 2: must be 0 or 1');
    end
else
    z1 = [zcrp(x); zcrp(-x)];
end

% find crossings when x==0 exactly
z0 = find( (x(1:length(x))==0) & ([x(2:length(x));0] ~= 0));

% concatenate and sort the two types of zero-crossings
i = sort([z0; z1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zz = zcrp(xx)  %only used in zcr
% find positive-going zero-crossing
z1 = find(diff(sign(xx)) == -2);
% find which out of current sample or next sample is closer to zero
[m, z2] = min([abs(xx(z1)), abs(xx(z1+1))], [], 2);
zz =  z1 -1 + z2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [frob]=frobfun(sp,p,m,offset)

% [frob]=frobfun(sp,p,m)
% 
% sp is the speech signal assumed to be preemphasised
% p  is the prediction order  : recomended to be 1 ms in above paper
% m  is the window length     : recomended to be 1 ms in above paper
% offset is shift for assymetric speech shape at start of voiced cycle -
% default 1.5ms.
%
% This function implements the frobenius norm based measure C defined in [4] below.
% It equals the square of the Frobenius norm of the m by p+1 data matrix divided by p+1
%
% Reference:
%   [4]  C. Ma, Y. Kamp, and L. F. Willems, “A Frobenius norm approach to glottal closure detection
%        from the speech signal,” IEEE Trans. Speech Audio Processing, vol. 2, pp. 258–265, Apr. 1994.


%   Revision History: 
%       V1.0 July 12th 2002:
%            Nov. 6th 2002 : if statement added to remove "last" midpoint 

%force p m and offset to be integers
p=round(p);
m=round(m);
offset=round(offset);

w=(p+1)*ones(1,m+p);
w(1:p)=1:p;
w(m+1:p+m)=p:-1:1;

w=w./(p+1); 
frob=filter(w,1,sp.^2);
frob(1:(round((p+m-1)/2) + offset))=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function goi=simplegci2goi(gci,pr)

% Estimate glottal opening instants by assuming a fixed closed-phase fraction

gci=round(gci);
maxpitch=max(medfilt1(diff(gci),7));

% calculate opening instants
for kg=1:length(gci)-1
    goi(kg)=gci(kg)+min(pr*(gci(kg+1)-gci(kg)),pr*maxpitch);
end;
kg=kg+1;
goi(kg)=round(gci(kg)+pr*(gci(kg)-gci(kg-1)));  %use the previous pitch period instead
goi=round(goi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tew,sew,y,toff]=xewgrdel(u,fs)

% implement EW group delay epoch extraction

dy_gwlen=voicebox('dy_gwlen');          % group delay evaluation window length
dy_fwlen=voicebox('dy_fwlen');          % window length used to smooth group delay

% perform group delay calculation

gw=2*floor(dy_gwlen*fs/2)+1;            % force window length to be odd
ghw=window('hamming',gw,'s');
ghw = ghw(:);                           % force to be a column (dmb thinks window gives a row - and he should know as he wrote it!)
ghwn=ghw'.*(gw-1:-2:1-gw)/2;            % weighted window: zero in middle

u2=u.^2;
yn=filter(ghwn,1,u2);
yd=filter(ghw,1,u2);
yd(abs(yd)<eps)=10*eps;                 % prevent infinities
y=yn(gw:end)./yd(gw:end);               % delete filter startup transient
toff=(gw-1)/2;
fw=2*floor(dy_fwlen*fs/2)+1;            % force window length to be odd
if fw>1
    daw=window('hamming',fw,'s');
    y=filter(daw,1,y)/sum(daw);         % low pass filter 
    toff=toff-(fw-1)/2;
end
[tew,sew]=zerocros(y,'n');              % find zero crossings

tew=tew+toff;                           % compensate for filter delay and frame advance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cfn=fnrg(gcic,frob,fs)

%Frobenious Energy Cost

dy_fxminf=voicebox('dy_fxminf');
frob=frob(:)';
mm=round(fs/dy_fxminf);
mfrob=maxfilt(frob,1,mm);
mfrob=[mfrob(floor(mm/2)+1:end) max(frob(end-ceil(mm/2):end))*ones(1,floor(mm/2))];
rfr=frob./mfrob;
Cfn=0.5-rfr(round(gcic));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gci=dpgci(gcic, s, Ch, fnwav, fs)

%DPGCI   Choose the best Glottal Closure Instances with Dynamic Programming
%   gci=dpgci(gcic, s(:), Ch, fnwav, fs) returns vectors of sample indices corresponding
%   to the instants of glottal closure in the speech signal s at sampling frequency fs Hz.
%
%   Inputs:
%   gcic    is a matrix whos first column are the glottal closure instance candidates and
%           the second column is 1 if the corresponding gci is derived from a zero crossing 
%           but zero if the gci is from a a projected zero crossing
%   s       is the speech signal - MUST be a column vector
%   Ch      the phase slope cost of every candidate
%   fnwav   is the frobenious norm function of s
%   fs      is the sampling frequncy
%
%   Outputs:
%   gci     is a vector of glottal closure instances chosen by the DP



%   Revision History: 
%   Bugs:  Constants are hardwired but defined in a structure like pv (defined in grpdelpv)
%         

% get algorithm parameters from voicebox()

dy_fxmin=voicebox('dy_fxmin');        % min larynx frequency (Hz)
dy_fxmax=voicebox('dy_fxmax');        % min larynx frequency (Hz)
dy_xwlen=voicebox('dy_xwlen');        % cross-correlation length for waveform similarity (sec)
dy_nbest=voicebox('dy_nbest');        % Number of NBest paths to keep
dy_spitch=voicebox('dy_spitch');              % scale factor for pitch deviation cost
wproj=voicebox('dy_cproj');           % cost of projected candidate
dy_cspurt=voicebox('dy_cspurt');           % cost of a talkspurt
dy_wpitch=voicebox('dy_wpitch');           % DP pitch weighting
dy_wener=voicebox('dy_wener');           % DP energy weighting
dy_wslope=voicebox('dy_wslope');           % DP group delay slope weighting
dy_wxcorr=voicebox('dy_wxcorr');           % DP cross correlation weighting



%Constants
Ncand=length(gcic);
sv2i=-(2*dy_spitch^2)^(-1);                 % scale factor for pitch deviation cost
nxc=ceil(dy_xwlen*fs);                      % cross correlation window length in samples
% === should delete any gci's that are within this of the end.
% === and for energy window

%Limit the search:
qrmin=ceil(fs/dy_fxmax);
qrmax=floor(fs/dy_fxmin);

%Cost and tracking r = current, q = previous, p = preprevious
cost=zeros(Ncand, dy_nbest); cost(:,:)=inf;    %Cost matrix, one row for each candidate
maxcost=zeros(Ncand,1); maxcost(:,:)=inf;   %Maximum cost in each row
imaxcost=ones(Ncand,1);                     %Index of maximum cost

prev = ones(Ncand, dy_nbest);                  %index of previous, q candidates
ind = ones(Ncand, dy_nbest);                   %index of p in row q (from prev)
qbest = [zeros(Ncand,1), ones(Ncand,2)]; % the minimum cost in any previous q [cost,q,i]

Cfn=fnrg(gcic(:,1),fnwav,fs);  %Frob.Energy Cost

%Add start and end state
% === should probably delete candidates that are too close to either end of the input
% === why do we ever need the additional one at the tail end ?
gcic=[[gcic(1,1)-qrmax-2 0];gcic;[gcic(end,1)+qrmax+2 0]];
Cfn=[0 Cfn 0];
Ch = [0 Ch 0];

% first do parallelized version


% rather complicated window specification is for compatibility with DYPSA 2
% === +1 below is for compatibility - probably a bug
wavix=(-floor(nxc/2):floor(nxc/2)+1)';                 % indexes for segments [nx2,1]
nx2=length(wavix);
sqnx2=sqrt(nx2);
g_cr=dy_wener*Cfn+dy_wslope*Ch+wproj*(1-gcic(:,2))';  % fixed costs

g_n=gcic(:,1)';                  % gci sample number [1,Ncand+2]
g_pr=gcic(:,2)';                 % unprojected flag [1,Ncand+2]
g_sqm=zeros(1,Ncand+1);         % stores: sqrt(nx2) * mean for speech similarity waveform
g_sd=zeros(1,Ncand+1);         % stores: 1/(Std deviation * sqrt(nx2)) for speech similarity waveform
f_pq=zeros((Ncand+1)*dy_nbest,1);   % (q-p) period for each node
f_c=repmat(Inf,(Ncand+1)*dy_nbest,1);    % cumulative cost for each node - initialise to inf
f_c(1)=0;                       % initial cost of zero for starting node
% f_costs=zeros(Ncand*dy_nbest,6);   % === debugging only remember costs of candidate
f_f=ones((Ncand+1)*dy_nbest,1);    % previous node in path
f_fb=ones((Ncand+1),1);    % points back to best end-of-spurt node
fbestc=0;                       % cost of best end-of-spurt node

qmin=2;
for r=2:Ncand+1   
%     if r==86
%         r;
%     end
    r_n=g_n(r);             % sample number of r = current candidate
    rix=dy_nbest*(r-1)+(1:dy_nbest);    % index range within node variables
    
    % determine the range of feasible q candidates
    qmin0=qmin;
    qmin=find(g_n(qmin0-1:r-1)<r_n-qrmax);      % qmin is the nearest candidate that is >qrmax away
    qmin=qmin(end)+qmin0-1;             % convert to absolute index of first viable candidate
    qmax=find(g_n(qmin-1:r-1)<=r_n-qrmin);      % qmax is the nearest candidate that is >=qrmin away
    qmax=qmax(end)+qmin-2;
    
    
    % calculate waveform similarity cost measure statistics
    
    sr=s(r_n+wavix);        % note s MUST be a column vector so sr is also
    wsum=sum(sr);
    g_sqm(r)=wsum/sqnx2;                % mean * sqrt(nx2)
    g_sd(r)=1/sqrt(sr.'*sr-wsum^2/nx2);   % 1/(Std deviation * sqrt(nx2))
    
    % now process the candidates
    
    if qmin<=qmax
        qix=qmin:qmax;      % q index
        nq=length(qix);
        % === should integrate the -0.5 into dy_wxcorr
        % === the factor (nx2-1)/(nx2-2) is to compensate for a bug in swsc()
        q_cas=-0.5*(nx2-1)/(nx2-2)*dy_wxcorr*(sum(s(repmat(g_n(qix),nx2,1)+repmat(wavix,1,nq)).*repmat(sr,1,nq),1)-g_sqm(qix)*g_sqm(r)).*g_sd(qix)*g_sd(r);
        % compare: i=35; Ca=swsc(g_n(qix(i)),g_n(r),s,fs); [i qix(i) r  g_n(qix(i)) g_n(r) dy_wxcorr*Ca q_cas(i)]
        
        % now calculate pitch deviation cost
        
        fix = 1+(qmin-1)*dy_nbest:qmax*dy_nbest;    % node index range
        f_qr=repmat(r_n-g_n(qix),dy_nbest,1);    % (r-p) period for each node
        f_pr=f_qr(:)+f_pq(fix);
        % === could absorb the 2 into sv2i
        f_nx=2-2*f_pr./(f_pr+abs(f_qr(:)-f_pq(fix)));
        f_cp=dy_wpitch*(0.5-exp(sv2i*f_nx.^2));
        % === fudge to match dypsa2.4 - could more efficiently be added
        % === onto the cost of a talkspurt end
        % === should be a voicebox parameter anyway
        f_cp(f_pq(fix)==0)=dy_cspurt*dy_wpitch;
        
        % now find the N-best paths
        
        [r_cnb,nbix]=sort(f_c(fix)+f_cp+reshape(repmat(q_cas,dy_nbest,1),nq*dy_nbest,1));
        f_c(rix)=r_cnb(1:dy_nbest)+g_cr(r);     % costs
        f_f(rix)=nbix(1:dy_nbest)+(qmin-1)*dy_nbest;       % traceback nodes
        f_pq(rix)=f_qr(nbix(1:dy_nbest));       % previous period
        % === f_costs is only for debugging
%         r;
%         f_costs(rix,1)=f_c(fix(nbix(1:dy_nbest)));
%         f_costs(rix,2)=wproj*(1-gcic(r,2));
%         f_costs(rix,3)=f_cp(nbix(1:dy_nbest));
%         f_costs(rix,4)=dy_wener*Cfn(r);
%         f_costs(rix,5)=dy_wslope*Ch(r);
%         f_costs(rix,6)=reshape(q_cas(1+floor((nbix(1:dy_nbest)-1)/dy_nbest)),dy_nbest,1);
        
        % check cost of using this candidate as the start of a new spurt
        % ==== the qmin>2 condition is for compatibility with dypsa 2 and
        % prevents any spurts starting until at least qrmax past the first
        % gci. This is probably a bug (see again below)
        iNb=rix(end);        
        if (qmin>2) && (f_c(f_fb(qmin-1))+wproj*(1-gcic(r,2))<f_c(iNb))        % compare with worst of Nbest paths
            f_f(iNb)=f_fb(qmin-1);
            % === for now we exclude the energy and phase-slope costs for compatibility with dypsa2
            % === this is probably a bug
            f_c(iNb)=f_c(f_fb(qmin-1))+wproj*(1-gcic(r,2));     % replace worst of the costs with start voicespurt cost
            f_pq(iNb)=0;                    % false pq period
        end
        if f_c(rix(1))<fbestc
            f_fb(r)=rix(1);                          % points to the node with lowest end-of-spurt cost
            % === should compensate for the pitch period cost incurred at the start of the next spurt
            % === note that a node can never be a one-node voicespurt on its own unless dy_nbest=1
            % since the start voices[purt option replaced the worst Nbest cost. This is probably good but
            % is a bit inconsistent.
            fbestc=f_c(rix(1));
        else
            f_fb(r)=f_fb(r-1);
        end
    else            % no viable candidates - must be the start of a voicespurt if anything
        % === for now we exclude the energy and phase-slope costs for compatibility with dypsa2
        % === this is probably a bug
        % ==== the qmin>2 condition is for compatibility with dypsa 2 and
        % prevents any spurts starting until at least qrmax past the first
        % gci. This is probably a bug (see again above)
        if (qmin>2)
            f_c(rix(1))=f_c(f_fb(qmin-1))+wproj*(1-gcic(r,2));  % cost of new voicespurt
            f_f(rix)=f_fb(qmin-1);                              % traceback to previous talkspurt end
            f_pq(rix)=0;                                        % previous period
        end
        f_fb(r)=f_fb(r-1);                                  % cannot be the end of a voicespurt
    end
end

% now do the traceback

gci = zeros(1,Ncand+1);

% === for compatibility with dypsa2, we force the penultimate candidate to be accepted
% === should be: i=f_fb(Ncand+1) but instead we pick the best of the penultimate candidate
i=rix(1)-dy_nbest;
if f_c(i-dy_nbest+1)<f_c(i)     % check if start of a talkspurt
    i=i-dy_nbest+1;
end
k=1;
while i>1
    j=1+floor((i-1)/dy_nbest);          % convert node number to candidate number
    gci(k)=g_n(j);
    i=f_f(i);
    k=k+1;
end
gci=gci(k-1:-1:1);           % put into ascending order 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


