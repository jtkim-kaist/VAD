function [xx,ii,m,v,mr,vr]=psycest(iq,x,r,xp,lf)
% Estimate multiple psychometric functions
%
% Usage: [xx,ii,m,v]=psycest(-n,p,q,xp,lf) % initialize n models
%        [xx,ii,m,v]=psycest(i,x,r)     % supply a trial result to psycest
%        [xx,ii,m,v]=psycest(i,x,r,o)   % supply a trial result to psycest and plot
%                    psycest(i,o)       % plot pdf of model i with plot options o
%        [xx,ii,m,v,mr,vr]=psycest(i,o) % Determine robust outputs in addition expcluding outliers [takes longer]
%              [p,q,msr]=psychest(0)    % output model parameters (or print them if no outputs)
%
% Inputs (see usage examples for argument positions):
%         -n        minus the number of models
%          p(:,n)   parameters for each model
%                      1  thresh [0.75]
%                      2  miss [0.04]
%                      3  guess [0.1]
%                      4  SNR min [-20]
%                      5  SNR max [20]
%                      6  Slope min (but over-ridden if <q.sl) [0]
%                      7  slope max [0.5]
%          q(:)     parameters common to all models (vector or struct)
%                      1  q.nx  number of SNR values in pdf [40]
%                      2  q.ns  number of slope values in pdf [21]
%                      3  q.nh  number of probe SNR values to evaluate [30]
%                      4  q.cs  weighting of slope relative to SRT in cost function [1]
%                      5  q.dh  minimum step size in dB for probe SNRs [0.2]
%                      6  q.sl  min slope at threshold (must be >0) [0.005]
%                      7  q.kp  number of std deviations of the pdfs to keep [4]
%                      8  q.hg  amount to grow expected gains in ni trials [1.3]
%                      9  q.cf  cost function: 1=variance, 2=entropy [2]
%                     10  q.pm  psychometric model: 1=logistic, 2=cumulative gaussian [1]
%                     11  q.lg  use log slope in pdf: 0=no, 1=yes [1]
%                     12  q.pp  Number of prior standard deviations in initial semi-range [1]
%                     13  q.pf  Probability floor (integrated over entire grid) [0.0001]
%                     14  q.ts  Number of std devs to explore [2]
%                     15  q.dp  Maximum probe SNR shift (dB) [10]
%                     16  q.it  Grid interpolation threshold [0.5]
%                     17  q.at  Axis change threshold [0.1]
%                     18  q.la  Look 2-ahead when choosing probe [1]
%                     19  q.op  Outlier probability [0.01]
%                     20  q.rx  Minimum range factor per iteration [0.5]
%          xp{n}(:) list of available probe SNR values
%          lf       log file ID
%          i        model probed
%          x        probe SNR value used
%          r        response: 0=wrong, 1=right.
%          o        plot option string:
%                   'p' Plot pdf
%                   'h' Plot probe history
%                   'x' Plot expected cost function for probe SNRs
%                   'c' Plot cost function evolution
%
% Outputs:
%          xx       recommended probe SNR
%          ii       recommended model to probe next
%          m(2,n,3) estimated srt and slope of all models
%                   m(:,:,1:3) are respectively the mean, mode (MAP) and marginal mode estimates
%          v(3,n)   estimated covariance matrix entries:
%                   [var(srt) cov(srt,slope) var(slope)]'
%          mr(2,n,3)robust estimated srt and slope of all models
%                   m(:,:,1:3) are respectively the mean, mode (MAP) and marginal mode estimates
%          vr(3,n)  robust estimated covariance matrix entries:
%                   [var(srt) cov(srt,slope) var(slope)]'
%          msr(:,3)  List probe snrs and results: [model probe-snr result]
%
% Algorithm parameters:
%
% The algorithm estimates the psychometric function for a number of models simultaneously.
% A psychometric function specifies the probability of correct recognition as a function of
% SNR and is completely specified by two free parameters: the SNR (in dB) and the slope (in 1/dB)
% at a specified target recognition probability (e.g. 0.5 or 0.75). The p(:,n) parameters specify
% some of the details of the psychometric function and can be different for each model:
%   p(1,n) gives the target recognition probability
%   p(2,n) gives the probability of incorrect recognition at very good SNRs (the "miss" probability).
%          If this value is made too low, a single unwarrented recognition error will have a large
%          effect of the estimated parameters and greatly lengthen the time to convergence. The default
%          value is 4%.
%   p(3,n) gives the probability of correct recognition at very poor SNRs (the "guess" probabiity).
%          This should be set to 1/N where N is the number of possible responses to a stimulus.
% p(4:5,n) These give the initial min and max SNR at the target recognition probability. They will
%           be adjusted by the algorithm if necessary but wrong values will delay convergence.
% p(6:7,n) These given the initial min and max slope (in probability per dB) at the target
%          recognition probability.
% The remaining parameters are shared between all the models and control how the algorithm operates.
%   q(1:2) These parameters determine the sampling size of the joint pdf in SNR and Slope respectively.
%   q(3)   This parameter specifies how many potential values of probe SNR to evaluate in order to
%          determine which is likely to give the most improvement in the parameter estimates.
%   q(4)   This parameter gives the relative weight of SNR and Slope in the cost function. Increasing
%          its value will improve the slope estimate (or log(slope) estimate) at the expense of the
%          SNR estimate. Actually its value is not all that critical.
%   q(5)   At each iteration psycest evaluates several potential probe SNR levels to see which is
%          likely to give the most improvement to the parameter estimates. This parameter fixes the
%          minimum spacing between these potential probe values. It only has an effect when the variances
%          of the parameter estimates are very small.
%   q(6)   To prevent the routine getting lost, this parameter specifies the smallest reasonable value
%          of the Slope parameter at the target recognition probability. Under normal circumstances, it
%          has no effect.
%   q(7)   For each model, the routine maintains a sampled version of the joint pdf of the SNR and slope.
%          The sampling grid is reclculated after each iteration to encompass this number of standard
%          deviations in each dimension.
%   q(8)   At each iteration, psycest advises which model to probe next on the basis of which
%          gives the greatest expected reduction in cost function. To ensure that all models
%          are periodically sampled, this expected reduction of each unprobed model is multiplied
%          by this parameter at each iteration. Thus a model will eventually be probed even if
%          its expected cost factor improvement is small.
%   q(9)   This determines whether the recommended probe SNR to use at the next iteration is chosen to
%          minimize the expected variance or the entropy of the estimated SNR and slope distributions.
%          My experience is that entropy (the default) gives faster convergence.
%  q(10)   This selects whether the underlying model is a logistic function (1) or a cumulative
%          gaussian (2). The choice makes little difference unless the "miss" or "guess" probabilities
%          are very very small.
%  q(11)   This selects whether the internal pdf samples "slope" (0) or "log(slope)" (1). It is recommended
%          that you stick to the default of log(slope) since this results in more symmetrical distributions. [1]
%  q(12)   This sets the number of std deviations of the prior corresponding to the limits of the
%          initial distribution set in p(4:7,n). A small value corresponds to a very weak prior. [1]
%  q(13)   Probability floor (integrated over entire grid) [0.0001]
%  q(14)   Number of std devs to explore when evaluating possible probe SNRs [2]
%  q(15)   Maximum probe SNR shift beyond the range of previous probes (dB) [10]
%  q(16)   Grid interpolation threshold. The grid is only interpolated if
%          it has shrunk by less that this factor [0.5]
%  q(17)   If the desired grid limits change by less than this fraction of
%          the grid length, they will be left unaltered. [0.1]
%  q(18)   Look 2-ahead when choosing probe [1]
%  q(19)   Outlier probability; probe results with lower probability than
%          this will be ignored [0.01]
%  q(20)   The ranges of SRT and log slope will not be reduced per iteration below this factor [0.5]

%      Copyright (C) Mike Brookes and Clement Doire 2009-2016
%      Version: $Id: psycest.m 9215 2016-12-17 22:37:49Z dmb $
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bugs/Suggestions:
% (2)use a structure for input parameters
% (3) entropy seems to change signifiacntly when the grid is rescaled
% (4)add a forgetting factor for old measurements
% (8) use quadratic interpolation to find cost function minimum
% (10) optionally output the whole pdf + its axis values
% (11) optionally output all model probe values
% (13) Should perhaps estimate and return the mean and slope compensated for the guess rate
%      i.e. if you want p, you actually test at guess+ p*(1-guess-miss)/(1-miss)
% (17) remember probe snrs, models and results and recalculate the entire
%      array when changing precision
% (18) could force probe SNRs to be a multiple of minimum step size
% (19) use a non-uniform prior e.e. 1/(1+x^2)
% (20) possibly use different parametrization (e.g. 1/slope or a second SRT threshold)
% (22) save inputs so that pdfs can be recalculated when rescaling
% (24) Parametrize slope as x=(100s^2-1)/20s to make the distribution more uniform
%      inverse is s=(sqrt(x^2+1)-x)/10; alternatively use log(slope)
% (25) optionally have no prior (i.e. maximum likelihood)
% (26) Check why the estimated variances are too small
% (27) Adapt range of potential probes if optimum is at the limit
% (28) Resample the pdf on the basis of the marginal distributions
% (29) Select the probe range on the basis of marginal distributions
% (30) Resample the pdfs by a factor of 2 always
% (31) use uneven samples in pdf concentrated in the middle
% (32) Use a parametric mixture including one-sided constants for the pdf
% (33) Selection of subset of prescribed probe SNRs is not optimum
% (34) use quadratic interpolation to select the probe SNR unless using fixed values
% (35) expand the pdf grid based on the effective number of samples (like particle filters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent variables
% sizes: nx=SNR values in pdf
%        ns=slope values in pdf (1/std dev of cumulative gaussian)
%        nh=potential test SNRs
%        ni=simultaneous models being estimated
% wq(ns*nx,ni)  = log prob of model; each column is a vectorized ns*nx matrix
% nr(1,10)      = parameter values: [SNR-pdf slope-pdf SNR-probe ...]
% pr(7,ni)      = input model parameters
% qr(4,ni)      = derived model parameters
% xq(nr(1),ni)  = SNR values in pdf
% sq(nr(2),ni)  = slope values in pdf
% mq(2,ni,3)    = estimated srt and slope of all models
% vq(3,ni)      = estimated covariance matrix entries:[var(srt) cov(srt,slope) var(slope)]'
% xn(1,ni)      = next probe value to use
% hn(1,ni)      = expected decrease in cost function after next probe
% hfact
% xz
% res(nres,7)   = results [model-number, probe-SNR, result]
% nres          = total number of probe results
% nresq(1,ni)   = number of probe results for each model
% xmm(2,ni)     = [min(probe-SNR); max(probe-SNR)]
% mq0(2,ni)      = prior means [x-mean; s-mean]
% pq0(2,ni)      = scale factors for prior distribution
% wfl           = floor relative to peak of log-pdf array
% sqmin         = minimum value of slope or log-slope
% LOG           = file ID of logging file
% mqr(2,ni,3)   = robust mean estimates (excluding outliers)
% vqr(3,ni)     = robust variance estimates (excluding outliers)
% nresr         = nres value at which robust estimates were last calculated

persistent wq xq sq nr pr qr mq vq xn hn hfact xz res nres nresq xmm mq0 pq0 wfl sqmin LOG mqr vqr nresr xlim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Initialization (iq<0)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iq<0  % initialization
    if nargin>4 && ~isempty(lf)
        LOG=lf;
        fprintf(LOG,'******************************\npsycest Initialization: %s\n******************************\n',datestr(now));
    else
        LOG=0;
    end
    ni=-iq;                 % number of models
    nres=0;                 % total number of probe-results so far
    nresq=zeros(1,ni);      % number of probe-results for each model
    res=zeros(1,7);         % array for saving results [expands as needed]
    pr=repmat([0.75 0.04 0.1 -20 20 0 0.5]',1,ni);          % default parameters
    if nargin>1
        if size(x,2)>1
            if size(x,2)>ni
                error('initialization parameter argument has too many columns');
            end
            pr(1:size(x,1),1:size(x,2))=x;
        else
            pr(1:size(x,1),:)=repmat(x,1,ni);
        end
    end
    nr=[40 21 30 1 0.2 0.02 4 1.3 2 1 1 1 0.0001 2 10 0.5 0.1 1 0.01 0.5]';                      % default parameter values
    nrf={'nx','ns','nh','cs','dh','sl','kp','hg','cf','pm','lg','pp', 'pf', 'ts', 'dp', 'it','at','la','op','rx'};     % parameter field names
    numnr=length(nr);
    if nargin>2 && ~isempty(r) % replace defaults with any parameter values specified in the call
        if isstruct(r)
            fnn=fieldnames(r);
            for i=1:length(fnn)
                mk=strcmp(fnn{i},nrf);
                if any(mk)
                    nr(mk)=r.(fnn{i});
                end
            end
        else
            nr(1:min(numnr,numel(r)))=r(:); % if parameters are specified as a vector, copy them across
        end
        nr(1:3)=round(nr(1:3));             % first three parameters must be integers
    end
    pr(6,:)=max(pr(6,:),nr(6));            	% low limit of slope in prob/dB
    nxq=nr(1);
    nsq=nr(2);
    nsxq=nxq*nsq;
    xq=(0:nxq-1)'*(pr(5,:)-pr(4,:))/(nxq-1)+repmat(pr(4,:),nxq,1);
    if nr(11)  % if log slope
        sqmin=log(nr(6)); % low limit for axis
        sq=(1-nsq:0)'*(log(pr(7,:))-log(max(pr(6,:),nr(6))))/(nsq-1)+repmat(log(pr(7,:)),nsq,1);
    else
        sqmin=nr(6); % low limit for axis
        sq=(1-nsq:0)'*(pr(7,:)-max(pr(6,:),nr(6)))/(nsq-1)+repmat(pr(7,:),nsq,1);
    end
    mq0=[mean(xq,1); mean(sq,1)];     % find mean of prior
    pq0=-2*nr(12)^2*[xq(end,:)-xq(1,:); sq(end,:)-sq(1,:)].^(-2); % scale factor for prior distribution
    wq=repmat(pq0(2,1)*(sq(:,1)-mq0(2,1)).^2,1,nxq)+repmat(pq0(1,1)*(xq(:,1)-mq0(1,1))'.^2,nsq,1); % log probs of prior (same for all models)
    wfl=log(nr(13)/nsxq); % floor relative to peak of wq array
    wq=repmat(wq(:),1,ni);          	% initialize to +-q.pp std deviations of gaussian prior
    qr=zeros(5,ni);
    qr(1,:)=1-pr(2,:)-pr(3,:);         	% prob range covered by cumulative gaussian
    if any(qr(1,:)<=0)
        [dum,i]=min(qr(1,:));
        error('Model %d: guess prob (%.2f) + miss prob (%.2f) is >=1',i,pr(3,i),pr(2,:));
    end
    qr(2,:)=(pr(1,:)-pr(3,:))./qr(1,:);                 % cumulative gaussian probability at threshold
    if any(abs(qr(2,:)-0.5)>=0.5)
        [dum,i]=max(qr(2,:)-0.5);
        error('Model %d: target SRT threshold (%.2f) must lie in the range guess prob (%.2f) to (1-miss prob) (%.2f)',i,pr(1,i),pr(3,i),1-pr(2,:));
    end
    switch nr(10)                                       % switch on model type
        case 1                                          % logistic model
            qr(3,:)=log(qr(2,:)./(1-qr(2,:)));        	% x position of target in std measure
            qr(4,:)=qr(1,:).*qr(2,:).*(1-qr(2,:));    	% slope*"stddev" at threshold
        case 2                                          % cumulative gaussian model
            qr(3,:)=norminv(qr(2,:));                   % x position of target in std measure
            qr(4,:)=qr(1,:).*normpdf(qr(3,:));       	% slope*stddev at threshold
        otherwise
            error('Unrecognised psychometric model selection');
    end
    if nr(9)< 1 || nr(9)>2
        error('Unrecognised cost function option');
    end
    mq=repmat(mq0,[1,1,3]);        % initial means, joint mode and marginal mode all are equal
    vq=[var(xq,1,1); zeros(1,ni); var(sq,1,1)];         % initial variances (actually ignores prior probs)
    mqr=mq;         % robust means and modes
    vqr=vq;         % robust variances
    nresr=0;        % robust calculation time
    xlim=[-Inf(1,ni); Inf(1,ni)]; % SNR limits for outliers
    %     hn=[1-nr(4) 0 nr(4)]*vq;  % artificially high value of cost function ensures all models are probed early on
    hn=Inf(1,ni);           % very high initial cost function
    hfact=nr(8)^(1/ni);     % growth factor to ensure that no model is neglected for too long
    xn=mq0(1,:);             % start at mean value
    xmm=[xn; xn];  % min/max of probe values for each model
    if nargin>=4 && ~isempty(xp)
        if iscell(xp)
            xz=xp;
        else
            xz=repmat(num2cell(xp(:)',2),ni,1);        % else replicate for each model
        end
        for i=1:ni
            [dum,j]=min(abs(xz{i}-mq0(1,i)));      % select the closest available probe to the mean
            xn(i)=xz{i}(j);
        end
    else
        xz=cell(ni,1);          % make empty cells
    end
elseif iq>0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     iq>0 means update or plot model iq           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        first check the plotting options          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<3
        if nargin==2
            po=x; % plot options are 2nd argument
        elseif ~nargout
            po='p';
        else
            po='';
        end
    else
        if nargin>3
            po=xp; % plot options are 4th argument
        else
            po='';
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     now get parameters of current model (iq)     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nxq=nr(1);                  % number of SNR values
    nsq=nr(2);                  % number of slope values
    nxh=nr(3);                  % number of probe-SNR values
    nsxq=nxq*nsq;               % size of log pdf array
    thresh=pr(1,iq);            % target probability at threshold
    guess=pr(3,iq);             % guess rate (1/choices)
    pscale=qr(1,iq);            % prob range left after subtracting miss and guess probs
    xtstd=qr(3,iq);             % x position of target in std measure
    xqi=xq(:,iq);               % SRT values of the pdf array
    sqi=sq(:,iq);               % slope (or log slope) values in PDF
    wqi=wq(:,iq);              	% log pdf array
    mqi= mq(:,iq,1);            % [xe; se] means
    vqi=vq(:,iq);               % [xv; sxv; sv]  covariance matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     rescale the pdfs if necessary         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ssd=sqrt(vqi(3));           % std deviation of slope
    xsd=sqrt(vqi(1));           % std deviation of SRT
    % determine range of new grid
    xqrange=xqi(end)-xqi(1);    % range of old SRT grid
    sqrange=sqi(end)-sqi(1);    % range of old slope grid
    if  nresq(iq)<2              % keep old limits if nresq(iq)<2
        xq2lim=xqi([1 end]);
        sq2lim=sqi([1 end]);
    else
        xqsemirange=max(nr(7)*xsd,0.5*nr(20)*xqrange);
        xq2lim=[mqi(1)-xqsemirange, mqi(1)+xqsemirange];
        sqsemirange=max(nr(7)*ssd,0.5*nr(20)*sqrange);
        sq2lim=[max(sqmin,mqi(2)-sqsemirange),mqi(2)+sqsemirange];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % To avoid unnecessary changes, keep old limits if the change  %
        % is less than nr(17) times the previous range [nr(17)=0.1]    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if abs(xq2lim(1)-xqi(1))<nr(17)*xqrange
            xq2lim(1)=xqi(1);
        end
        if abs(xq2lim(2)-xqi(end))<nr(17)*xqrange
            xq2lim(2)=xqi(end);
        end
        
        if abs(sq2lim(1)-sqi(1))<nr(17)*sqrange
            sq2lim(1)=sqi(1);
        end
        if abs(sq2lim(2)-sqi(end))<nr(17)*sqrange
            sq2lim(2)=sqi(end);
        end
    end
    xq2=linspace(xq2lim(1),xq2lim(2),nxq)';   % new x axis values
    sq2=linspace(sq2lim(1),sq2lim(2),nsq)';   % new s axis values
    wqup=2;                                   % update flag
    if xq2(1)<xqi(1) || xq2(end)>xqi(end) || sq2(1)<sqi(1) || sq2(end)>sqi(end)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if extrapolating, recalculate log-pdfs from saved data  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if LOG, fprintf(LOG,'N=%d:%d, extrapolate: x(%.2g,%.2g)->(%.2g,%.2g) and s(%.2g,%.2g)->(%.2g,%.2g)\n',iq,nresq(iq),xqi([1 end]),xq2([1 end]),sqi([1 end]),sq2([1 end])); end
        wq2=repmat(pq0(2,iq)*(sq2-mq0(2,iq)).^2,1,nxq)+repmat(pq0(1,iq)*(xq2-mq0(1,iq))'.^2,nsq,1); % log probs of prior for model iq
        ires=find(res(1:nres,1)==iq); % find the results for this model
        if nr(11)                   % if log slope
            sqis=exp(sq2)/qr(4,iq);	% inverse std dev of gaussian (proportional to slope)
        else
            sqis=sq2/qr(4,iq);     	% inverse std dev of gaussian (proportional to slope)
        end
        switch nr(10)   % switch on model type
            case 1
                for i=1:length(ires)
                    j=ires(i);
                    r0=res(j,3)==0;
                    wq2=wq2+log(r0+(1-2*r0)*(guess+pscale*((1+exp(sqis*xq2'-xtstd-repmat(sqis,1,nxq)*res(j,2))).^(-1)))); %  P(l | r,x)
                end
            case 2
                for i=1:length(ires)
                    j=ires(i);
                    r0=res(j,3)==0;
                    wq2=wq2+log(r0+(1-2*r0)*(guess+pscale*normcdf(repmat(sqis,1,nxq)*res(j,2)-sqis*xq2'+xtstd))); %  P(l | r,x)
                end
        end
    else                                                    % possibly do quadratic interpolation
        wq2=max(reshape(wqi,nsq,nxq)-max(wqi(:)),wfl);    	% turn into a normalized, clipped matrix for easy interpolation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % use quadratic interpolation in SRT axis %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (xq2(end)-xq2(1))/(xqrange)>nr(16)      	% if range has shrunk by < nr(16), leave the SRT axis unchanged
            xq2=xqi;                              	% copy the old SRT axis
            wqup=1;                                   % update flag
        else                                                % else do quadratic interpolation
            if LOG, fprintf(LOG,'N=%d:%d, interpolate x: x(%.2g,%.2g)->(%.2g,%.2g)\n',iq,nresq(iq),xqi([1 end]),xq2([1 end])); end
            xqf=(xq2-xqi(1))/(xqi(2)-xqi(1));
            xqj=ceil(xqf);
            xqf=xqj-xqf;
            xqg=1-xqf;
            xqh=0.25*xqf.*xqg; % Quadratic coefficient
            xqf((xqj<=0) | (xqj>nxq))=0;
            xqg((xqj<0) | (xqj>=nxq))=0;
            xqh((xqj<1) | (xqj>=nxq-1))=0;
            wq2=wq2(:,min(max(xqj,1),nxq)).*repmat(xqf'+xqh',nsq,1)+wq2(:,min(max(xqj+1,1),nxq)).*repmat(xqg'+xqh',nsq,1) - ...
                (wq2(:,min(max(xqj-1,1),nxq))+wq2(:,min(max(xqj+2,1),nxq))).*repmat(xqh',nsq,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % use quadratic interpolation in slope axis %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (sq2(end)-sq2(1))/(sqrange)>nr(16)       % if range has shrunk by < nr(16), leave the slope axis unchanged
            sq2=sqi;                                % copy the old slope axis
            wqup=wqup-1;                            % update flag
        else
            if LOG, fprintf(LOG,'N=%d:%d, interpolate s: s(%.2g,%.2g)->(%.2g,%.2g)\n',iq,nresq(iq),sqi([1 end]),sq2([1 end])); end
            sqf=(sq2-sqi(1))/(sqi(2)-sqi(1));
            sqj=ceil(sqf);
            sqf=sqj-sqf;
            sqg=1-sqf;
            sqh=0.25*sqf.*sqg; % Quadratic coefficient
            sqf((sqj<=0) | (sqj>nsq))=0;
            sqg((sqj<0) | (sqj>=nsq))=0;
            sqh((sqj<1) | (sqj>=nsq-1))=0;
            wq2=wq2(min(max(sqj,1),nsq),:).*repmat(sqf+sqh,1,nxq)+wq2(min(max(sqj+1,1),nsq),:).*repmat(sqg+sqh,1,nxq) - ...
                (wq2(min(max(sqj-1,1),nsq),:)+wq2(min(max(sqj+2,1),nsq),:)).*repmat(sqh,1,nxq);
        end
    end
    if wqup>0                               % check if we have updated either axis
        wq2=max(wq2(:)-max(wq2(:)),wfl);    % turn back into a normalized, clipped vector
        sqi=sq2;                            % update slope (or log slope) values in PDF
        xqi=xq2;                            % update SRT values of the pdf array
        wqi=wq2;                            % log pdf array
        sq(:,iq)=sqi;                       % save new s axis (log slope or slope)
        xq(:,iq)=xqi;                       % save new x axis (SRT)
        wq(:,iq)=wqi;                       % save log-pdf
    end
    if nr(11)                   % if log slope
        sqis=exp(sqi)/qr(4,iq);	% inverse std dev of gaussian (proportional to slope)
    else
        sqis=sqi/qr(4,iq);     	% inverse std dev of gaussian (proportional to slope)
    end
    if nargin>=3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update pdfs with a new probe result  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nres>=size(res,1)        % ensure there is space for a new result
            res=[res; zeros(nres,7)]; % double the size of the results array
        end
        nres=nres+1;                % increment the total number of results
        nresq(iq)=nresq(iq)+1;        % increment the number of results for this model
        res(nres,1:3)=[iq,x,r];       % save the latest result
        xmm(:,iq)=[min(x,xmm(1,iq)); max(x,xmm(2,iq))]; % update min/max probe values
        hn=hn*hfact;         % increase entropy gains to ensure all models get a chance
        
        %     disp(res(1:nres,:));  % ********** diagnostic
        
        %
        % update log probabilities with the previous test result
        % If r==1, we multiply (add in log) by a horizontally reflected version of the psychometric function that equals
        % the threshold (e.g. 0.5 or 0.75) at the probe SNR, x, for all slopes.
        % If r==0, we multiply (add in log) by 1- this reflected version which therefore equals (1-thresh) at x.
        %
        r0=r==0;
        switch nr(10)   % switch on model type
            case 1
                wqi=wqi+log(r0+(1-2*r0)*(guess+pscale*((1+exp(reshape(sqis*xqi'-xtstd,nsxq,1)-repmat(sqis,nxq,1)*x)).^(-1)))); %  P(l | r,x)
            case 2
                wqi=wqi+log(r0+(1-2*r0)*(guess+pscale*normcdf(repmat(sqis,nxq,1)*x-reshape(sqis*xqi'-xtstd,nsxq,1)))); %  P(l | r,x)
            otherwise
                error('Unrecognised psychometric model selection');
        end
        wq(:,iq)=wqi;             % save updated probabilities
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate mean and covariance and entropy  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ewqi=exp(wqi-max(wqi));             % unnormalized probability vector
    ewqi=ewqi/sum(ewqi);
    wqsx=reshape(ewqi,nsq,nxq);         % normalized probabilities
    px=sum(wqsx,1);                             %  p(x0)
    ps=sum(wqsx,2);                             %  p(s0)
    xe=px*xqi;                                  % E(x0)
    se=ps'*sqi;                                 % E(s0)
    
    [pxpk,xm]=max(px);                          % marginal mode in x
    if xm>1 && xm<nxq                           % use quadratic interpolation in log prob if possible
        [dum,xm2]=quadpeak(log(px(xm-1:xm+1))');
        xm=xm+xm2-2;
    end
    xm=(2-xm)*xqi(1)+(xm-1)*xqi(2);             % marginal mode(x0)
    
    [pspk,sm]=max(ps);
    if sm>1 && sm<nsq                           % use quadratic interpolation in log prob if possible
        [dum,sm2]=quadpeak(log(ps(sm-1:sm+1)));
        sm=sm+sm2-2;
    end
    sm=(2-sm)*sqi(1)+(sm-1)*sqi(2);             % marginal mode(s0)
    
    [wqpk,j]=max(ewqi);
    i=1+floor((j-1)/nsq);
    j=j-nsq*(i-1);
    if i>1 && i<nxq && j>1 && j<nsq             % use quadratic interpolation in log prob if possible
        [dum,ji]=quadpeak(wqi(repmat((j-1:j+1)',1,3)+repmat(nsq*(i-2:i),3,1)));
        i=i+ji(2)-2;
        j=j+ji(1)-2;
    end
    xj=(2-i)*xqi(1)+(i-1)*xqi(2);               % joint mode  x
    sj=(2-j)*sqi(1)+(j-1)*sqi(2);               % joint mode: s
    xv=px*(xqi.^2)-xe^2;                        % Var(x0)
    sv=ps'*(sqi.^2)-se^2;                       % Var(s0)
    sxv=ewqi'*(repmat(sqi-se,nxq,1).*reshape(repmat(xqi'-xe,nsq,1),nsxq,1)); % Cov(s0*x0)
    mq(:,iq,1)=[xe; se];                        % save means
    mq(:,iq,2)=[xj; sj];                        % save joint mode
    mq(:,iq,3)=[xm; sm];                        % save marginal modes
    vq(:,iq)=[xv; sxv; sv];                     % save covariance matrix
    xh=(px*log(px)')*(xqi(1)-xqi(2));           % differential entropy h(x0)
    sh=(ps'*log(ps))*(sqi(1)-sqi(2));          	% differential entropy h(s0)
    if nargin==3
        switch nr(9) % cost function: 1=variance, 2=entropy
            case 1
                res(nres,4:7)=[xe se xv sv];           	% save info for plotting history
            case 2
                res(nres,4:7)=[xe se xh sh];                                       	% find the minimum of cost function
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now calculate robust mean and covariance for all models   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargout>4 || ~isempty(po)   	% need to update robust statistics if >4 outputs or if plotting
        if  nr(19)==0               % if no outliers allowed in any case
            mr=mq;
            vr=vq;
        elseif nresr==nres          % if robust statistics are all up-to-date
            mr=mqr;
            vr=vqr;
        else                        % update robust statistics for all models
            for jq=1:length(nresq)             % loop for each model
                if any(res(nresr+1:nres,1)==jq)   % if we have had any recent probe results for this model
                    guessj=pr(3,jq);                % guess rate (1/choices)
                    pscalej=qr(1,jq);            % prob range left after subtracting miss and guess probs
                    xtstdj=qr(3,jq);               % x position of target in std measure
                    xqj=xq(:,jq);               % SRT values of the pdf array
                    sqj=sq(:,jq);               % slope (or log slope) values in PDF
                    % first determine the SNR range to include
                    xej=mq(1,jq,1);  % mean SRT
                    sej=mq(2,jq,1);   % mean slope or log slope
                    resj=res(res(:,1)==jq,:);  % results for this model
                    if nr(11)                   % if log slope
                        sqisj=exp(sej)/qr(4,jq);	% inverse std dev of gaussian (proportional to slope) at mean slope
                    else
                        sqisj=sej/qr(4,jq);     	% inverse std dev of gaussian (proportional to slope) at mean slope
                    end
                    % find upper bound for SNR with r=0
                    switch nr(10)   % switch on model type
                        case 1
                            xlim(2,jq)=xej-(xtstdj+log(nr(19)/(pscalej-nr(19))))/sqisj; % outlier if r=0 and x > xlim
                        case 2
                            xlim(2,jq)=xej+(norminv((pscalej-nr(19))/pscalej)-xtstdj)/sqisj; % outlier if r=1 and x < xlim
                    end
                    mou=resj(:,2)>=xlim(2,jq) & resj(:,3)==0;       % find high outliers
                    if nr(19)>guessj  % find lower bound for SNR with r=1
                        switch nr(10)   % switch on model type
                            case 1
                                xlim(1,jq)=xej-(xtstdj+log(pscalej/(nr(19)-guessj)-1))/sqisj; % outlier if r=1 and x < xlim
                            case 2
                                xlim(1,jq)=xej+(norminv((nr(19)-guessj)/pscalej)-xtstdj)/sqisj; % outlier if r=1 and x < xlim
                        end
                        mou=mou | resj(:,2)<=xlim(1,jq) & resj(:,3)==1;       % add in low outliers
                    end
                    if any(mou) % if any results must be excluded
                        % now recalculate pdf excluding bad points
                        if nr(11)                   % if log slope
                            sqisj=exp(sqj)/qr(4,jq);	% inverse std dev of gaussian (proportional to slope) for all grid points
                        else
                            sqisj=sqj/qr(4,jq);     	% inverse std dev of gaussian (proportional to slope) for all grid points
                        end
                        wqj=repmat(pq0(2,jq)*(sqj-mq0(2,jq)).^2,1,nxq)+repmat(pq0(1,jq)*(xqj-mq0(1,jq))'.^2,nsq,1); % log probs of prior
                        switch nr(10)               % switch on model type
                            case 1                  % nr(10)=1: logistic model
                                for j=find(mou==0)'
                                    r0=resj(j,3)==0;
                                    wqj=wqj+log(r0+(1-2*r0)*(guessj+pscalej*((1+exp(sqisj*xqj'-xtstdj-repmat(sqisj,1,nxq)*resj(j,2))).^(-1)))); %  P(l | r,x)
                                end
                            case 2                  % nr(10)=2: cumulative gaussian model
                                for j=find(mou==0)'
                                    r0=resj(j,3)==0;
                                    wqj=wqj+log(r0+(1-2*r0)*(guessj+pscalej*normcdf(repmat(sqisj,1,nxq)*resj(j,2)-sqisj*xqj'+xtstdj))); %  P(l | r,x)
                                end
                        end
                        % now calculate robust means and modes
                        wqj=wqj(:);                         % turn 2D array into a column vector
                        ewqj=exp(wqj-max(wqj));             % unnormalized probability vector
                        ewqj=ewqj/sum(ewqj);
                        wqsxr=reshape(ewqj,nsq,nxq);         % normalized probabilities
                        pxr=sum(wqsxr,1);                             %  p(x0)
                        psr=sum(wqsxr,2);                             %  p(s0)
                        xer=pxr*xqj;                                  % E(x0)
                        ser=psr'*sqj;                                 % E(s0)
                        [pxpk,xmr]=max(pxr);                          % marginal mode in x
                        if xmr>1 && xmr<nxq                           % use quadratic interpolation in log prob if possible
                            [dum,xm2]=quadpeak(log(pxr(xmr-1:xmr+1))');
                            xmr=xmr+xm2-2;
                        end
                        xmr=(2-xmr)*xqj(1)+(xmr-1)*xqj(2);             % marginal mode(x0) (requires the xqi to be uniformly spaced)
                        [pspk,smr]=max(psr);
                        if smr>1 && smr<nsq                           % use quadratic interpolation in log prob if possible
                            [dum,sm2]=quadpeak(log(psr(smr-1:smr+1)));
                            smr=smr+sm2-2;
                        end
                        smr=(2-smr)*sqj(1)+(smr-1)*sqj(2);        	% marginal mode(s0)
                        [wqpk,j]=max(ewqj);                         % find max of 1-D array
                        i=1+floor((j-1)/nsq);                       % convert to 2-D index (i,j)
                        j=j-nsq*(i-1);
                        if i>1 && i<nxq && j>1 && j<nsq             % use quadratic interpolation in log prob if possible
                            [dum,ji]=quadpeak(wqj(repmat((j-1:j+1)',1,3)+repmat(nsq*(i-2:i),3,1)));
                            i=i+ji(2)-2;
                            j=j+ji(1)-2;
                        end
                        xjr=(2-i)*xqj(1)+(i-1)*xqj(2);           	% joint mode  x
                        sjr=(2-j)*sqj(1)+(j-1)*sqj(2);           	% joint mode: s
                        xvr=pxr*(xqj.^2)-xer^2;                   	% Var(x0)
                        svr=psr'*(sqj.^2)-ser^2;                  	% Var(s0)
                        sxvr=ewqj'*(repmat(sqj-ser,nxq,1).*reshape(repmat(xqj'-xer,nsq,1),nsxq,1)); % Cov(s0*x0)
                        mqr(:,jq,:)=reshape([xer xjr xmr; ser sjr smr],[2 1 3]); % save means, joint modes and marginal modes
                        vqr(:,jq)=[xvr; sxvr; svr];                	% save covariance matrix
                    else   % if there are no outliers, just copy non-robust results
                        mqr(:,jq,:)= mq(:,jq,:);                    % use non-robust means
                        vqr(:,jq)= vq(:,jq);                      	% and non-robust variances
                    end
                end
            end
            nresr=nres;                 % all robust means/modes are now up-to-date
            mr=mqr;                     % send to output arguments
            vr=vqr;
        end
        if nr(11)    % if using log-slope
            mr(2,:,:)=exp(mr(2,:,:));  % convert to real slope
            vr(2,:)=vr(2,:).*mr(2,:,1); % correct the covariance
            vr(3,:)=vr(3,:).*mr(2,:,1).^2; % and the slope variance
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now estimate the next probe SNR            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~numel(xz{iq})                               % if no list of probe SNRs was specified
        ytry=exp(0.25i*pi*(0:7))';                  % points around the circle
        ytry=[real(ytry) imag(ytry)];
        [vtry,dtry]=eig([xv sxv; sxv sv]);          % eigendecomposition of covariance matrix
        tryxs=repmat([xe,se],8,1)+nr(14)*ytry*sqrt(dtry)*vtry';
        pmin=0.05;                                  % target probe success probability
        if nr(11)
            tryxs(:,2)=qr(4,iq)*exp(-tryxs(:,2)); 	% convert log(slope) to std dev
        else
            tryxs(:,2)=qr(4,iq)*tryxs(:,2).^(-1); 	% convert slope to std dev
        end
        dp=nr(15);  % maximum shift of probe value outside previous range
        switch nr(10)   % switch by model type
            case 1 % logistic
                qmax=min(xmm(2,iq)+dp,max(tryxs(:,1)+(log((1-pmin)/pmin)-xtstd)*tryxs(:,2)));
                qmin=max(xmm(1,iq)-dp,min(tryxs(:,1)+(log(pmin/(1-pmin))-xtstd)*tryxs(:,2)));
            case 2 % cumulative gaussian
                qmax=min(xmm(2,iq)+dp, max(tryxs(:,1)+(norminv(1-pmin)-xtstd)*tryxs(:,2)));
                qmin=max(xmm(1,iq)-dp,min(tryxs(:,1)+(norminv(pmin)-xtstd)*tryxs(:,2)));
        end
        dxt=max(nr(5),(qmax-qmin)/nxh);             % minimum step size of nr(5) [0.2 dB]
        xt=(qmin+qmax)/2+((1:nxh)-(1+nxh)/2)*dxt;
    else                                            % if a specific list of probe SNRs exists
        xzi=xz{iq};                                 % xzi is the list of available probe SNRs
        if numel(xzi)<=nxh                          % use all available probe SNRs if there are not too many
            xt=xzi;
        else
            [xt,ixt]=min(abs(xzi-xe));                                  % find the closest one to xe ** not necessarily optimum ***
            ixt=max(1,min(1+numel(xzi)-nxh,ixt-floor((1+nxh)/2)));      % arrange symmetrically around xt
            xt=xzi(ixt:min(ixt+nxh-1,numel(xzi)));
        end
    end
    nxhp=length(xt);  % xt are the potential probe SNRs
    % Now find the probe value that minimizes the expected value of the cost function
    % In the following: l = parameters of psychometric function, x = the probe SNR and r = the probe result
    switch nr(10)
        case 1
            prt=guess+pscale*((1+exp(repmat(reshape(sqis*xqi'-xtstd,nsxq,1),1,nxhp)-repmat(sqis,nxq,1)*xt)).^(-1)); %  P(r=1 | l,x)
        case 2
            prt=guess+pscale*normcdf(repmat(sqis,nxq,1)*xt-repmat(reshape(sqis*xqi'-xtstd,nsxq,1),1,nxhp)); %  P(r=1 | l,x)
    end
    wqt=repmat(ewqi,1,nxhp);
    hminr=zeros(2,1);                   % space for minimum expected cost function for each r0
    hminj=zeros(1,nxhp);   % space for minimum expected cost function for each x0
    if nr(18) % if doing look 2-ahead
        pl1=prt.*wqt;                       % posterior prob given success = p(l | x0,r0=1) unnormalized
        pl0=wqt-pl1;                        % posterior prob given failure = p(l | x0,r0=0) unnormalized
        pr0=sum(pl1,1);                     % p(r0=1 | x0)=Sum{P(r0=1 | l,x0)*P(l)} [note each column of wqt is normalized] (row vector)
        plxr=reshape([pl0./repmat(1-pr0,nsxq,1) pl1./repmat(pr0,nsxq,1)],nsxq,nxhp,2); % posterior prob p(l | x0,r0) column-normalized
        nx0=nxhp; % outer loop for each x0
        nr0=2; % inner loop for each r0
    else  % if only doing look 1-ahead
        nx0=1;  % only execute outer loop once
        nr0=1;  % only execute inner loop once
        pr0=0; % dummy value (used at end of outer loop)
    end
    hx2=zeros(nx0,nxhp,nr0);        % space for square array of expected cost functions (for possible plotting only)
    for j=1:nx0                     % loop for each possible probe SNR, x0, (or only once is look-1-ahead)
        for jr=1:nr0             	% loop for each possible test result, r0=jr-1 (or only once is look-1-ahead)
            if nr(18) % if doing look 2-ahead
                wqt=repmat(plxr(:,j,jr),1,nxhp);         % posterior prob p(l | x0=xt(j),r0=jr-1) column-normalized
            end
            pl1=prt.*wqt;                     % posterior prob given success = p(l | x,r=1) unnormalized
            pl0=wqt-pl1;                      % posterior prob given failure = p(l | x,r=0) unnormalized
            prh=sum(pl1,1);                     % p(r | x)=Sum{P(r | l,x)*P(l)} [note each column of wqtjr is normalized] (row vector)
            pl1=pl1./repmat(prh,nsxq,1);        % posterior prob given success = p(l | x,r=1) normalized
            pl0=pl0./repmat(1-prh,nsxq,1);      % posterior prob given failure = p(l | x,r=0) normalized
            
            px1=squeeze(sum(reshape(pl1,nsq,nxq,[]),1));    % p(x0 | x,r=1)
            px0=squeeze(sum(reshape(pl0,nsq,nxq,[]),1));    % p(x0 | x,r=0)
            ps1=squeeze(sum(reshape(pl1,nsq,nxq,[]),2));    % p(s0 | x,r=1)
            ps0=squeeze(sum(reshape(pl0,nsq,nxq,[]),2));    % p(s0 | x,r=0)
            xet1=xqi'*px1;                                  % E(x0 | x,r=1)
            xvt1=(xqi.^2)'*px1-xet1.^2;                     % Var(x0 | x,r=1)
            xet0=xqi'*px0;                                  % E(x0 | x,r=0)
            xvt0=(xqi.^2)'*px0-xet0.^2;                     % Var(x0 | x,r=0)
            xvt=xvt1.*prh+xvt0.*(1-prh);                    % E(Var(x0 | x ))
            set1=sqi'*ps1;                                  % E(s0 | x,r=1)
            svt1=(sqi.^2)'*ps1-set1.^2;                     % Var(s0 | x,r=1)
            set0=sqi'*ps0;                                  % E(s0 | x,r=0)
            svt0=(sqi.^2)'*ps0-set0.^2;                     % Var(s0 | x,r=0)
            svt=svt1.*prh+svt0.*(1-prh);                    % E(Var(s0 | x ))
            xht1=sum(log(px1).*px1,1);                      % -H(x0 | x,r=1)
            xht0=sum(log(px0).*px0,1);                      % -H(x0 | x,r=0)
            xht=(xht1.*prh+xht0.*(1-prh))*(xqi(1)-xqi(2));	% h(x0 | x)
            sht1=sum(log(ps1).*ps1,1);                      % -H(s0 | x,r=1)
            sht0=sum(log(ps0).*ps0,1);                      % -H(s0 | x,r=0)
            sht=(sht1.*prh+sht0.*(1-prh))*(sqi(1)-sqi(2));	% h(s0 | x)
            switch nr(9) % cost function: 1=variance, 2=entropy
                case 1
                    hx=(xvt + nr(4)*svt)/(1+nr(4));           	% expected cost function for each possible test SNR
                case 2
                    hx=(xht + nr(4)*sht)/(1+nr(4));            	% expected cost function for each possible test SNR                                      	% find the minimum of cost function
            end
            hx2(j,:,jr)=hx; % save expected cost function for possible plotting
            [hminr(jr),ix]=min(hx);       % find the minimum of cost function for each value of r (0 or 1)
            %         fprintf('Probe range: %.3g %.3g; choose %.3g\n',xt(1),xt(end),xt(ix));
        end
        hminj(j)=[1-pr0(j) pr0(j)]*hminr; % expected cost function for each x0 assuming x1 is chosen optimally
    end
    if nr(18) % if doing look 2-ahead
        [hminr(1),ix]=min(hminj);                         % find the minimum of cost function
    end
    % now ix indexes the optimal probe choice and hminr(1) is the resultant expected cost function
    xn(iq)=xt(ix);  % save the next probe snr for this model
    % calculate the expected decrease in cost function (to decide which model to probe next)
    switch nr(9)  % cost function: 1=variance, 2=entropy
        case 1                                              % minimize variance
            hn(iq)=(xv + nr(4)*sv)/(1+nr(4))-hminr(1);    	% save the expected decrease in cost function for this model
        case 2                                              % minimize entropy
            hn(iq)=(xh + nr(4)*sh)/(1+nr(4))-hminr(1);    	% save the expected decrease in cost function for this model
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  now do plotting                                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(po)
        h=gcf;
        if isreal(h)
            fig0=round(h); % get figure number
        else
            fig0=get(h,'number');  % in new versions of matlab use this method
        end
        axm=[3 -2; 2 -1]; % matrix to extend an axis beackwards by two steps
        for i=1:length(po)
            poi=po(i);
            switch poi
                case 'p'    % plot pdf
                    nxq=nr(1);
                    nsq=nr(2);
                    xqj=[axm*xqi(1:2);xqi];
                    sqj=[axm*sqi(1:2);sqi];
                    figure(fig0);
                    imagesc(xqj,sqj,[zeros(1,2) px/pxpk; zeros(1,nxq+2);ps/pspk zeros(nsq,1) wqsx/wqpk]);
                    hold on
                    plot(xe,se,'ko',xm,sm,'k^',xj,sj,'k*');
                    plot(xqj(2),se,'wo',xqj(2),sm,'w^',xqj(2),sj,'w*');
                    plot(xe,sqj(2),'wo',xm,sqj(2),'w^',xj,sqj(2),'w*');
                    if any(mqr(:,iq,:)~=mq(:,iq,:))
                        plot(mqr(1,iq,1),mqr(2,iq,1),'go',mqr(1,iq,3),mqr(2,iq,3),'g^',xj,sj,'g*');
                        plot(xqj(end),mqr(2,iq,1),'go',xqj(end),mqr(2,iq,3),'g^',xqj(end),mqr(2,iq,2),'g*');
                        plot(mqr(1,iq,1),sqj(end),'go',mqr(1,iq,3),sqj(end),'g^',mqr(1,iq,2),sqj(end),'g*');
                        title('Joint pdf: o mean, * mode, \Delta marg mode (green=robust)');
                    else
                        title('Joint pdf: o mean, * mode, \Delta marginal mode');
                    end
                    t=linspace(0,2*pi,200);
                    xcir=cos(t);
                    scir=sin(t);
                    vcir=sqrt((sv*xcir.^2+xv*scir.^2-2*sxv*xcir.*scir)/(xv*sv-sxv^2));
                    plot(xe+xcir./vcir,se+scir./vcir,'w-');
                    hold off
                    axis 'xy';
                    %     colorbar;
                    %     cblabel('Relative Probability');
                    xlabel(sprintf('SNR @ %d%%SRT (dB)',round(pr(1)*100)));
                    if nr(11)
                        ylabel('Log psychometric Slope at threshold (prob/dB)');
                    else
                        ylabel('Psychometric Slope at threshold (prob/dB)');
                    end
                    fig0=fig0+1;
                case 'c' % plot cost evolution
                    resi=res(res(:,1)==iq,:);
                    nresi=size(resi,1);
                    if nresi>0
                        figure(fig0);
                        plot(1:nresi,resi(:,6)/(1+nr(4)),'--r',1:nresi,(resi(:,6) + nr(4)*resi(:,6))/(1+nr(4)),'-b');
                        axisenlarge([-1 -1.03]);
                        xlabel('Trial Number');
                        switch nr(9)  % cost function: 1=variance, 2=entropy
                            case 1                                              % minimize variance
                                ylabel('Weighted Variance: SRT and Total');
                            case 2                                              % minimize entropy
                                ylabel('Weighted Entropy: SRT and Total');
                        end
                        title ('Cost function evolution');
                        fig0=fig0+1;
                    end
                case 'h' % plot history
                    resi=res(res(:,1)==iq,:);
                    nresi=size(resi,1);
                    if nresi>0
                        if nr(11)
                            resi(:,5)=exp(resi(:,5)); % convert to real slope
                            ylim=xe+([-2 0 thresh 1 3]-thresh)/exp(se);
                        else
                            ylim=xe+([-2 0 thresh 1 3]-thresh)/se;
                        end
                        figure(fig0);
                        plot(1:nresi,resi(:,[4 4 4])+resi(:,5).^(-1)*([0 thresh 1]-thresh),'-b');
                        hold on;
                        % could plot outliers in a different colour by checking mou(res(:,1)==iq)
                        msky=resi(:,3)==1;
                        plot(find(msky),resi(msky,2),'r+');
                        msky=resi(:,3)==0;
                        plot(find(msky),resi(msky,2),'bo');
                        gcaylim=get(gca,'ylim');
                        poutb=0;
                        if nr(19)>0 && xlim(2,iq)<gcaylim(2) % plot upper outlier bound
                            plot([1;nresi],[1;1]*xlim(2,iq),'-r');
                            texthvc(nresi,xlim(2,iq),'max o','rtr');
                        end
                        if nr(19)>guess && xlim(1,iq)>gcaylim(1) % plot upper outlier bound
                            plot([1;nresi],[1;1]*xlim(1,iq),'-r');
                            texthvc(nresi,xlim(1,iq),'min+','rbr');
                        end
                        hold off;
                        xlabel('Trial Number');
                        ylabel('Probe SNR (dB) [o,x=result]');
                        axisenlarge(-1.02);
                        title('SRT mean + slope intercepts');
                        fig0=fig0+1;
                    end
                case 'x' % plot probe choice
                    figure(fig0);
                    if nr(18) % if doing look 2-ahead
                        xtaxm=xt(1:2)*axm';
                        hx2e=hx2(:,:,1).*repmat(1-pr0',1,nxhp)+hx2(:,:,2).*repmat(pr0',1,nxhp);
                        [dum,ix2m]=min(hx2,[],2); % ix2 indices to minimize hx2
                        % ==========  note that the following line assumes that the xt()
                        % are uniformly spaced (not necessarily true)
                        imagesc([xtaxm xt],xt,[hminj' repmat(max(hx2e(:)),nxhp,1) hx2e]);
                        axis xy;
                        colorbar;
                        hold on
                        [hnmin,ii]=max(hn);         % find the next model to probe
                        plot(xtaxm(1),xn(iq),'>w');
                        if iq==ii
                            plot(xt(1:2)*axm(1,:)',xn(iq),'*w'); % circle the next model to be probed
                        end
                        plot([xtaxm(1) xt(end)],xn(iq)*[1 1],':w');
                        plot(xt(ix2m(:,1,1)),xt,'ow');
                        plot(xt(ix2m(:,1,2)),xt,'+w');
                        hold off
                        xlabel('SNR probe 2 (dB) [+,o=best for R1=1,0]');
                        ylabel(sprintf('SNR probe 1 (dB) [> = %.2f]',xn(iq)));
                        title('Expected Cost after 2 probes');
                    else % if doing look 1-ahead only
                        plot(xt,hx2(1,:,1),'b');
                        axisenlarge([-1 -1.03]);
                        caylim=get(gca,'ylim');
                        hold on
                        plot(xn([iq iq]),caylim,':k');
                        xlabel(sprintf('Next probe SNR (dB) [min @ %.2f]',xn(iq)));
                        ylabel('Cost')
                        title('Expected Cost after next probe');
                    end
                    fig0=fig0+1;
            end
        end
    end
end
if iq~=0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  now select the appropriate model to probe next  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [hnmin,ii]=max(hn);         % chose model with the biggest expected decrease
    xx=xn(ii);
    m=mq;
    v=vq;
    if nr(11)    % if using log-slope
        m(2,:,:)=exp(m(2,:,:));  % convert to real slope
        v(2,:)=v(2,:).*m(2,:,1); % correct the covariance
        v(3,:)=v(3,:).*m(2,:,1).^2; % and the slope variance
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     iq=0 means output  model parameters and probe results   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xx=pr;
    ii=nr;
    m=res(1:nres,:);
    if ~nargout                     % print them if no output arguments to call
        fid=LOG+(LOG==0);           % use standard out if LOG==0
        pdesc={'Threshold','Miss prob','Guess Prob','Min SNR','Max SNR','Min Slope','Max Slope'};
        fprintf(LOG,'\n*********\nModel-specific Parameters\n');
        for i=1:7
            fprintf(LOG,'%4d) %s: ',i,pdesc{i});
            if size(pr,2)>1
                fprintf(LOG,'%.5g, ',pr(i,1:size(pr,2)-1));
            end
            fprintf(LOG,'%.5g\n',pr(i,size(pr,2)));
        end
        qdesc={'nx  SNR values in PDF', ...
            'ns  Slope values in PDF', ...
            'nh  Probe SNR values to evaluate', ...
            'cs  Weighting of slope relative to SRT in cost function', ...
            'dh  Min step size in dB for probe SNRs', ...
            'sl  Min slope at threshold', ...
            'kp  Std deviations of the pdfs to keep', ...
            'hg  Amount to grow expected gains in ni trials', ...
            'cf  Cost function', ...
            'pm  Psychometric model', ...
            'lg  Use log slope as parameter', ...
            'pm  Prior std devs in semi-width', ...
            'pf  Integrated probability floor', ...
            'ts  Number of std devs to explore', ...
            'dp  Max shift in probe value', ...
            'it  Grid interpolation threshold', ...
            'at  Grid boundary tolerance', ...
            'la  Look 2-ahead when choosing probe', ...
            'op  Oulier probability threshold', ...
            'rx  min grid shrink factor'};
        qoptf=[9 10]; % fields with options
        qoptval={'variance','entropy'; ...
            'logistic','cumulative gaussian'};
        fprintf(LOG,'\nShared Parameters\n');
        for i=1:length(nr)
            fprintf(LOG,'%4d) %s: ',i,qdesc{i});
            j=find(qoptf==i,1);
            if numel(j)
                fprintf(LOG,'%d=%s\n',nr(i),qoptval{j,nr(i)});
            else
                fprintf(LOG,'%.5g\n',nr(i));
            end
        end
        fprintf(LOG,'\n');
    end
end