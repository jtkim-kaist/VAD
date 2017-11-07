function [fx,tx,pv,fv]=fxpefac(s,fs,tinc,m,pp)
%FXPEFAC PEFAC pitch tracker [FX,TT,PV,FV]=(S,FS,TINC,M,PP)
%
% Input:   s(ns)      Speech signal
%          fs         Sample frequency (Hz)
%          tinc       Time increment between frames (s) [0.01]
%                     or [start increment end]
%          m          mode
%                     'g' plot graph showing waveform and pitch
%                     'G' plot spectrogram with superimposed pitch using
%                         options pp.sopt [default: 'ilcwpf']
%                     'x' use external files for algorithm parameter
%                         initialization: fxpefac_g and fxpefac_w
%          pp         structure containing algorithm parameters
%
% Outputs: fx(nframe)     Estimated pitch (Hz)
%          tx(nframe)     Time at the centre of each frame (seconds).
%          pv(nframe)     Probability of the frame of being voiced
%          fv             structure containing feature vectors
%                           fv.vuvfea(nframe,2) = voiced/unvoiced GMM features

% References
%  [1]  S. Gonzalez and M. Brookes. PEFAC - a pitch estimation algorithm robust to high levels of noise.
%       IEEE Trans. Audio, Speech, Language Processing, 22 (2): 518-530, Feb. 2014.
%       doi: 10.1109/TASLP.2013.2295918.
%  [2]  S.Gonzalez and M. Brookes,
%       A pitch estimation filter robust to high levels of noise (PEFAC), Proc EUSIPCO,Aug 2011.

% Bugs/Suggestions
% (1) do long files in chunks
% (2) option of n-best DP

%	   Copyright (C) Sira Gonzalez and Mike Brookes 2011
%      Version: $Id: fxpefac.m 5624 2015-01-12 12:15:34Z dmb $
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

persistent w_u m_u v_u w_v m_v v_v dpwtdef
% initialize persistent variables
if ~numel(w_u)

    % voiced/unvoiced decision based on 2-element feature vector
    % (a) mean power of the frame's log-freq spectrum (normalized so its short-term average is LTASS)
    % (b) sum of the power in the first three peaks
    %===== VUV
    if nargin>3 && any(m=='x')
        fxpefac_g;     % read in GMM parameters
        fxpefac_w;     % read in Weights parameters
    else
        w_u=[0.1461799 0.3269458 0.2632178 0.02331986 0.06360947 0.1767271 ]';

        m_u=[13.38533 0.4199435 ;
             12.23505 0.1496836 ;
             12.76646 0.2581733 ;
             13.69822 0.6893078 ;
             9.804372 0.02786567 ;
             11.03848 0.07711229 ];

        v_u=reshape([0.4575519 0.002619074 0.002619074 0.01262138 ;
             0.7547719 0.008568089 0.008568089 0.001933864 ;
             0.5770533 0.003561592 0.003561592 0.00527957 ;
             0.3576287 0.01388739 0.01388739 0.04742106 ;
             0.9049906 0.01033191 0.01033191 0.0001887114 ;
             0.637969 0.009936445 0.009936445 0.0007082946 ]',[2 2 6]);

        w_v=[0.1391365 0.221577 0.2214025 0.1375109 0.1995124 0.08086066 ]';

        m_v=[15.36667 0.8961554 ;
             13.52718 0.4809653 ;
             13.95531 0.8901121 ;
             14.56318 0.6767258 ;
             14.59449 1.190709 ;
             13.11096 0.2861982 ];

        v_v=reshape([0.196497 -0.002605404 -0.002605404 0.05495016 ;
             0.6054919 0.007776652 0.007776652 0.01899244 ;
             0.5944617 0.0485788 0.0485788 0.03511229 ;
             0.3871268 0.0292966 0.0292966 0.02046839 ;
             0.3377683 0.02839657 0.02839657 0.04756354 ;
             1.00439 0.03595795 0.03595795 0.006737475 ]',[2 2 6]);
    end
    %===== PDP
    %     dfm = -0.4238; % df mean
    %     dfv = 3.8968; % df variance (although treated as std dev here)
    %     delta = 0.15;
    %     dflpso=[dfm 0.5/(log(10)*dfv^2) -log(2*delta/(dfv*sqrt(2*pi)))/log(10)]; % scale factor & offset for df pdf
    %     dpwtdef=[1.0000, 0.8250, 1.3064, 1.9863]; % default DP weights
    dpwtdef=[1.0000, 0.8250, 0.01868, 0.006773, 98.9, -0.4238]; % default DP weights
    %===== END

end
% Algorithm parameter defaults

p.fstep=5;              % frequency resolution of initial spectrogram (Hz)
p.fmax=4000;            % maximum frequency of initial spectrogram (Hz)
p.fres = 20;            % bandwidth of initial spectrogram (Hz)
p.fbanklo = 10;         % low frequency limit of log filterbank (Hz)
p.mpsmooth = 21;       % width of smoothing filter for mean power
% p.maxtranf = 1000;      % maximum value of tranf cost term
p.shortut = 7;          % max utterance length to average power of entire utterance
p.pefact = 1.8;         % shape factor in PEFAC filter
p.numopt = 3;           % number of possible frequencies per frame
p.flim = [60 400];      % range of feasible fundamental frequencies (Hz)
p.w = dpwtdef;          % DP weights
% p.rampk = 1.1;          % constant for relative-amplitude cost term
% p.rampcz = 100;         % relative amplitude cost for missing peak
p.tmf = 2;              % median frequency smoothing interval (s)
p.tinc = 0.01;          % default frame increment (s)
p.sopt = 'ilcwpf';      % spectrogram options

% update parameters from pp argument

if nargin>=5 && isstruct(pp)
    fnq=fieldnames(pp);
    for i=1:length(fnq)
        if isfield(p,fnq{i})
            p.(fnq{i})=pp.(fnq{i});
        end
    end
end

% Sort out input arguments
if nargin>=3  && numel(tinc)>0
    p.tinc = tinc;   % 0.01 s between consecutive time frames
end
if nargin<4
    m='';
end

% Spectrogram of the mixture
fmin = 0; fstep = p.fstep; fmax = p.fmax;
fres = p.fres;  % Frequency resolution (Hz)
[tx,f,MIX]=spgrambw(s,fs,fres,[fmin fstep fmax],[],p.tinc);
nframes=length(tx);
txinc=tx(2)-tx(1);  % actual frame increment
%  ==== we could combine spgrambw and filtbankm into a single call to spgrambw or use fft directly ====
% Log-frequency scale
[trans,cf]=filtbankm(length(f),2*length(f)-1,2*f(end),p.fbanklo,f(end),'usl');
O = MIX*trans'; % Original spectrum in Log-frequency scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude Compression

% Calculate alpha based on LTASS ratios
ltass = stdspectrum(6,'p',cf);
auxf = [cf(1),(cf(1:end-1)+cf(2:end))./2,cf(end)];
ltass = ltass.*diff(auxf);                  % weight by bin width

% estimated ltass
O = O.*repmat(diff(auxf),nframes,1);     % weight spectrum by bin width
O1 = O;

if tx(end)<p.shortut                        % if it is a short utterance
    eltass = mean(O,1);                     % mean power per each frequency band
    eltass = smooth(eltass,p.mpsmooth);     % smooth in log frequency
    eltass= eltass(:).';                    % force a row vector

    % Linear AC
    alpha = (ltass)./(eltass);
    alpha = alpha(:).';
    alpha = repmat(alpha,nframes,1);
    O = O.*alpha;                           % force O to have an average LTASS spectrum

    % ==== should perhaps exclude the silent portions ***
else                                        % long utterance

    tsmo = 3; % time smoothing over 3 sec
    stt = round(tsmo/txinc);
    eltass = timesm(O,stt);
    eltass = smooth(eltass,p.mpsmooth);     % filter in time and log frequency

    % Linear AC
    alpha = repmat(ltass,nframes,1)./(eltass);
    O = O.*alpha;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the filter to detect the harmonics
ini = find(cf>3*cf(1));
sca = cf/cf(ini(1)); % bin frequencies start at approximately 0.33 with sca(ini(1))=1 exactly

% Middle
sca = sca(sca<10.5 & sca>0.5);  % restrict to 0.5 - 10.5 times fundamental

sca1 = sca;
filh = 1./(p.pefact-cos(2*pi*sca1));
filh = filh - sum(filh(1:end).*diff([sca1(1),(sca1(1:end-1)+sca1(2:end))./2,sca1(end)]))/sum(diff([sca1(1),(sca1(1:end-1)+sca1(2:end))./2,sca1(end)]));

posit = find(sca>=1);  % ==== this should just equal ini(1) ====
negat = find(sca<1);
numz = length(posit)-1-length(negat);
filh = filh./max(filh);
filh = [zeros(1,numz) filh]; % length is always odd with central value = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter the log-frequency scaled spectrogram
B = imfilter(O,filh);  % does a convolution with zero lag at centre of filh

% Feasible frequency range
numopt = p.numopt; % Number of possible fundamental frequencies per frame
flim = p.flim;
pfreq = find(cf>flim(1) & cf<flim(2)); % flim = permitted fx range = [60 400]
ff = zeros(nframes,numopt);
amp = zeros(nframes,numopt);
for i=1:nframes
    [pos,peak]=v_findpeaks(B(i,pfreq),[],5/(cf(pfreq(2))-cf(pfreq(1)))); % min separation = 5Hz @ fx=flim(1) (could pre-calculate) ====
    if numel(pos)
        [peak,ind]=sort(peak,'descend');
        pos = pos(ind);                     % indices of peaks in the B array
        posff = cf(pfreq(pos));             % frequencies of peaks
        fin = min(numopt,length(posff));
        ff(i,1:fin)=posff(1:fin);           % save both frequency and amplitudes
        amp(i,1:fin)=peak(1:fin);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probabilitly of the frame of being voiced

% voiced/unvoiced decision based on 2-element feature vector
% (a) mean power of the frame's log-freq spectrum (normalized so its short-term average is LTASS)
% (b) sum of the power in the first three peaks

pow = mean(O,2);

vuvfea = [log(pow) 1e-3*sum(amp,2)./(pow+1.75*1e5)];

% %%%%%%%%%%%%%%%%%%%%%

pru=gaussmixp(vuvfea,m_u,v_u,w_u);  % Log probability of being unvoiced
prv=gaussmixp(vuvfea,m_v,v_v,w_v);  % Log probability of being voiced

pv=(1+exp(pru-prv)).^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic programming

% w(1): relative amp, voiced local cost
% w(2): median pitch deviation cost
% w(3): df cost weight
% w(4): max df cost
% w(5): relative amp cost for missing peaks (very high)
% w(6): df mean

w = p.w;

% Relative amplitude 
camp = -amp./repmat(max(amp,[],2),1,numopt);  % relative amplitude used as cost
camp(amp==0)=w(5); % If no frequency found

% Time interval for the median frequency
tmf = p.tmf; % in sec
inmf = round(tmf/txinc);

%--------------------------------------------------------------------------
% FORWARDS
% Initialize values
cost = zeros(nframes,numopt);
prev = zeros(nframes,numopt);
medfx = zeros(nframes,1);
dffact=2/txinc;

% First time frame
% cost(1,:) = w(1)*ramp(1,:);
cost(1,:) = w(1)*camp(1,:);  % only one cost term for first frame
fpos = ff(1:min(inmf,end),1);
mf=median(fpos(pv(1:min(inmf,end))>0.6));   % calculate median frequency of first 2 seconds
if isnan(mf)
    mf=median(fpos(pv(1:min(inmf,end))>0.5));
    if isnan(mf)
        mf=median(fpos(pv(1:min(inmf,end))>0.4));
        if isnan(mf)
            mf=median(fpos(pv(1:min(inmf,end))>0.3)); % ==== clumsy way of ensuring that we take the best frames ====
            if isnan(mf)
                mf=0;
            end
        end
    end
end
medfx(1)=mf;

for i=2:nframes              % main dynamic programming loop
    if i>inmf
        fpos = ff(i-inmf:i,1);  % fpos is the highest peak in each frame
        mf=median(fpos(pv(1:inmf)>0.6));  % find median frequency over past 2 seconds
        if isnan(mf)
            mf=median(fpos(pv(1:inmf)>0.5));
            if isnan(mf)
                mf=median(fpos(pv(1:inmf)>0.4));
                if isnan(mf)
                    mf=median(fpos(pv(1:inmf)>0.3));% ==== clumsy way of ensuring that we take the best frames ====
                    if isnan(mf)
                        mf=0;
                    end
                end
            end
        end
    end
    medfx(i)=mf;
    % Frequency difference between candidates and cost
    df = dffact*(repmat(ff(i,:).',1,numopt) - repmat(ff(i-1,:),numopt,1))./(repmat(ff(i,:).',1,numopt) + repmat(ff(i-1,:),numopt,1));
    costdf=w(3)*min((df-w(6)).^2,w(4));

    % Cost related to the median pitch
    if mf==0                                   % this test was inverted in the original version
        costf = zeros(1,numopt);
    else
        costf = abs(ff(i,:) - mf)./mf;
    end
    [cost(i,:),prev(i,:)]=min(costdf + repmat(cost(i-1,:),numopt,1),[],2); % ==== should we allow the possibility of skipping frames ? ====
    cost(i,:)=cost(i,:)+w(2)*costf + w(1)*camp(i,:);  % add on costs that are independent of previous path

end

% Traceback

fx=zeros(nframes,1);
ax=zeros(nframes,1);
best = zeros(nframes,1);

nose=find(cost(end,:)==min(cost(end,:))); % ==== bad method (dangerous) ===
best(end)=nose(1);
fx(end)=ff(end,best(end));
ax(end)=amp(end,best(end));
for i=nframes:-1:2
    best(i-1)=prev(i,best(i));
    fx(i-1)=ff(i-1,best(i-1));
    ax(i-1)=amp(i-1,best(i-1));
end

if nargout>=4
    fv.vuvfea=vuvfea;  % voiced-unvoiced features
    fv.best=best;  % selected path
    fv.ff=ff;  % pitch candidates
    fv.amp=amp;  % pitch candidate amplitudes
    fv.medfx=medfx;  % median pitch
    fv.w=w;  % DP weights
    fv.dffact=dffact;  % df scale factor
    fv.hist = [log(mean(O,2)) sum(amp,2)./((mean(O,2)))];
end

if ~nargout || any(m=='g') || any(m=='G')
    nax=0;  % number of axes sets to link
    msk=pv>0.5; % find voiced frames as a mask
    fxg=fx;
    fxg(~msk)=NaN; % allow only good frames
    fxb=fx;
    fxb(msk)=NaN; % allow only bad frames
    if any(m=='G') || ~nargout && ~any(m=='g')
        clf;
        spgrambw(s,fs,p.sopt); % draw spectrogram with log axes
        hold on
        plot(tx,log10(fxg),'-b',tx,log10(fxb),'-r'); % fx track
        yy=get(gca,'ylim');
        plot(tx,yy(1)+yy*[-1;1]*(0.02+0.05*pv),'-k'); % P(V) track
        hold off
        nax=nax+1;
        axh(nax)=gca;
        if any(m=='g')
            figure;   % need a new figure if plotting two graphs
        end
    end
    if any(m=='g')
        ns=length(s);
        [tsr,ix]=sort([(1:ns)/fs 0.5*(tx(1:end-1)+tx(2:end))']); % intermingle speech and frame boundaries
        jx(ix)=1:length(ix); % create inverse index
        sp2fr=jx(1:ns)-(0:ns-1);  % speech sample to frame number
        spmsk=msk(sp2fr);   % speech sample voiced mask
        sg=s;
        sg(~spmsk)=NaN;   % good speech samples only
        sb=s;
        sb(spmsk)=NaN;    % bad speech samples only
        clf;
        subplot(5,1,1);
        plot(tx,pv,'-b',(1:ns)/fs,0.5*mod(cumsum(fx(sp2fr)/fs),1)-0.6,'-b');
        nax=nax+1;
        axh(nax)=gca;
        ylabel('\phi(t), P(V)');
        set(gca,'ylim',[-0.65 1.05]);
        subplot(5,1,2:3);
        plot((1:ns)/fs,sg,'-b',(1:ns)/fs,sb,'-r');
        nax=nax+1;
        axh(nax)=gca;
        subplot(5,1,4:5);
        plot(tx,fxg,'-b',tx,fxb,'-r');
        ylabel('Pitch (Hz)');
        %         semilogy(tx,fxg,'-b',tx,fxb,'-r');
        %         ylabel(['Pitch (' yticksi 'Hz)']);
        set(gca,'ylim',[min(fxg)-30 max(fxg)+30]);
        nax=nax+1;
        axh(nax)=gca;
    end
    if nax>1
        linkaxes(axh,'x');
    end
end

function y=smooth(x,n)
nx=size(x,2);
nf=size(x,1);
c=cumsum(x,2);
y=[c(:,1:2:n)./repmat(1:2:n,nf,1) (c(:,n+1:end)-c(:,1:end-n))/n (repmat(c(:,end),1,floor(n/2))-c(:,end-n+2:2:end-1))./repmat(n-2:-2:1,nf,1)];

function y=timesm(x,n)
if ~mod(n,2)
    n = n+1;
end
nx=size(x,2);
nf=size(x,1);
c=cumsum(x,1);
mid = round(n/2);
y=[c(mid:n,:)./repmat((mid:n).',1,nx); ...
    (c(n+1:end,:)-c(1:end-n,:))/n; ...
    (repmat(c(end,:),mid-1,1) - c(end-n+1:end-mid,:))./repmat((n-1:-1:mid).',1,nx)];
