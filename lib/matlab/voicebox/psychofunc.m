function p=psychofunc(m,q,x,r)
% Calculate psychometric functions: trial success probability versus SNR
%
% Usage: p=psychofunc('',q,x)       % calculate probabilities
%        b=psychofunc('r',q,x)       % generate boolean variables with success prob p
%        p=psychofunc(m,q,x,r)     % Calculate likelihoods for observations r
%        x=psychofunc([m 'i'],q,p) % Calculate inverse
%
% Inputs:
%          m        mode string [may be omitted if not required]
%                      'n'   do not normalize likelihoods
%                      'f'   do not squeeze output arrays to remove singleton dimensions
%                      'i'   calculate inverse function
%                      'r'   calculate binary random variables with probability p
%                      ['s'   calculate sweet points for threshold and slope]
%                      ['d'   calculate partial derivatives with respect to q(1:5)]
%                      'g'   plot graph
%                      'G'   plot image
%                      'c'   include colourbar
%          q        model parameters. Either a column vector with a single model,
%                   a matrix with one model per column or a cell array with multiple values for
%                   some or all of the parameters
%                      1  probability at threshold [0.5]
%                      2  threshhold [0 dB]
%                      3  slope at threshold [0.1 prob/dB ]
%                      4  miss or lapse probability [0]
%                      5  guess probability   [0]
%                      6  psychometric function type [1]
%                          1 = logistic
%                          2 = cumulative Gaussian
%                          3 = Weibull
%                          [4 = reversed Weibull]
%                          [5 = Gumbell]
%                          [6 = reversed Gumbell]
%          x        vector of SNR values
%          r        test results (0 or 1) corresponding to x
%          p        vector of probabilities
%
% Outputs:
%          p        array of probabilities or random variates ('r' option).
%                   p is a squeezed 7-dimensional array
%                   whose dimensions correspond to x followed by the 6 model parameter entries.
%                   if q is a cell array, singleton dimensions are removed unless the 'f' option is given.
%          x        Inverse function gives SNR, x, as a function of p
%          b        array of boolean variables

%      Copyright (C) Mike Brookes 2009-2010
%      Version: $Id: psychofunc.m 9302 2017-01-18 16:19:20Z dmb $
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

% first sort out input arguments
minp=0.01;          % minimum probability to use for inverse function by default
qq=[0.5 0 0.1 0 0 1]';  % default values for q
if nargin<4
    r=[];
    if nargin<3
        x=[];
        if nargin<2
            q=[];
            if ~nargin
                m='';
            end
        end
    end
end
if ~ischar(m);      % mode argument is optional
    r=x;
    x=q;
    q=m;
    m='';
end
sq=size(q);
ckmod=0;
if iscell(q)
    nq=ones(1,6);
    qax=num2cell([0; qq]);  % used for plotting
    for i=1:min(numel(q),6)
        nq(i)=numel(q{i});
        if nq(i)>=1
            nr=size(qq,2);
            qax{i+1}=q{i};
            if i<=5             % do not replicate for multiple models
                qq=repmat(qq,1,nq(i));
                qq(i,:)=reshape(repmat(q{i}(:)',nr,1),1,nr*nq(i));
            else
                qq(i,:)=repmat(q{i}(1),1,nr);
            end
        end
    end
    nq=max(nq,1);
    nmod=nq(6);
    if nmod>1      % list of models to use
        modlist=q{6};
    else
        modlist=qq(6,1);  % default model
    end
else
    nq=sq(2);
    if nq
        ql=repmat(qq,1,nq);
        ql(1:sq(1),:)=q;
    else
        ql=qq;
        nq=1;
    end
    modlist=unique(ql(6,:));
    nmod=length(modlist);
    ckmod=nmod>1;   % need to check model list
    qq=ql;
end
% now perform the calculation
nx=numel(x);
npt=50;       % number of points
if any(m=='i') % doing inverse
    if ~nx
        nx=npt;
        xlim=[max(qq(5,:)),1-max(qq(4,:))]*[1-minp minp; minp 1-minp];
        x=linspace(xlim(1),xlim(2),nx)';
    end
    p=zeros([nx nq]);  % space for SNRs
    ia=0;
    for i=1:nmod % loop for each model type
        mod=modlist(i);
        if ckmod
            qq=ql(:,ql(6,:)==mod);
        end
        pscale=1-qq(4,:)-qq(5,:);
        pstd=(qq(1,:)-qq(5,:))./pscale; % prob target compensating for miss and lapse probs
        sstd=qq(3,:)./pscale; % slope compensating for miss and lapse probs
        px=x(:)*pscale.^(-1)-repmat(qq(5,:)./pscale,nx,1);  % adjust for miss and lapse probs
        switch mod
            case 1
                beta=sstd./(pstd.*(1-pstd));
%                 alpha=qq(2,:)+log((1-pstd)./pstd)./beta;
                px=repmat(qq(2,:)+log((1-pstd)./pstd)./beta,nx,1)-log(px.^(-1)-1).*repmat(beta.^(-1),nx,1);
            case 2   % cumulative Gaussian function
                xtstd=norminv(pstd); % x position of target in std measure
                sig=normpdf(xtstd)./sstd;
                px= repmat(qq(2,:)-sig.*xtstd,nx,1) + repmat(sig,nx,1).*norminv(px);
            case 3
                wlog=log(1-pstd);
                kbeta=sstd./((pstd-1).*wlog);
                alpha=qq(2,:)-log(-wlog)./kbeta;
                px=repmat(alpha,nx,1)+log(-log(1-px)).*repmat(kbeta.^(-1),nx,1);
            otherwise
                error('Invalid psychometric model index');
        end
        if ckmod
            p(:,ql(6,:)==i)=px;
        else
            ib=ia+numel(p)/nmod;
            p(ia+1:ib)=px(:);
            ia=ib;
        end
    end
else % doing forward mapping
    if ~nx
        ef=2;         % expansion factor
        nx=npt;
        x=linspace(min(qq(2,:)-ef*(qq(1,:)-qq(5,:))./qq(3,:)), ...
            max(qq(2,:)+ef*(1-qq(1,:)-qq(4,:))./qq(3,:)),nx)';
    end
    p=zeros([nx nq]);  % space for probabilities
    ia=0;
    for i=1:nmod % loop for each model type
        mod=modlist(i);
        if ckmod
            qq=ql(:,ql(6,:)==mod);
        end
        pscale=1-qq(4,:)-qq(5,:);  % prob range excluding miss and lapse probs
        pstd=(qq(1,:)-qq(5,:))./pscale; % prob target compensating for miss and lapse probs
        sstd=qq(3,:)./pscale; % slope compensating for miss and lapse probs
        switch mod
            case 1   % logistic function
                beta=sstd./(pstd.*(1-pstd));
%                 alpha=qq(2,:)+log((1-pstd)./pstd)./beta;
                px=(1+exp(repmat(beta.*qq(2,:)+log((1-pstd)./pstd),nx,1)-x(:)*beta)).^(-1);
            case 2   % cumulative Gaussian function
                xtstd=norminv(pstd); % x position of target in std measure
                sigi=sstd./normpdf(xtstd);
                px=normcdf(x(:)*sigi-repmat(qq(2,:).*sigi-xtstd,nx,1));
            case 3
                wlog=log(1-pstd);
                kbeta=sstd./((pstd-1).*wlog);
                alpha=qq(2,:)-log(-wlog)./kbeta;
                px=1-exp(-exp(x(:)*kbeta-repmat(alpha.*kbeta,nx,1)));
            otherwise
                error('Invalid psychometric model index');
        end
        px=repmat(qq(5,:),nx,1)+repmat(pscale,nx,1).*px;  % adjust for miss and lapse probs
        if ckmod
            p(:,ql(6,:)==i)=px;
        else
            ib=ia+numel(p)/nmod;
            p(ia+1:ib)=px(:);
            ia=ib;
        end
    end
    if numel(r)                 % we are calculating likelihoods
        mk=r(:)==0;
        p(mk,:)=1-p(mk,:);      % invert probability for results that are zero
        if nx>1
            if any(m=='n')
                p=prod(p,1);
            else
                p=sum(log(p),1);
                p=exp(p-max(p(:)));
                p=p/sum(p(:));     % normalize to equal 1
            end
            nx=1;
        end
    end

end
pg=p;       % save unsqueezed p for plotting
if ~any(m=='f') && iscell(q) % remove all singleton dimensions
    szp=size(p);
    szq=szp(szp>1);
    szq=[szq ones(1,max(0,2-numel(szq)))];
    p=reshape(p,szq);
end
if any(m=='r') && ~any(m=='i');
    p=rand(size(p))<p;
end

if ~nargout || any(lower(m)=='g')
    clf;
    szp=[nx nq];
    czp=sum(szp>1);
    if czp>0  % check if there is anything to plot
        if iscell(q)
            axlab={'Input SNR','Threshold prob','Threshold SNR','Threshold slope','Lapse prob','Guess prob','Sigmoid type'};
            [szs,izs]=sort(szp,'descend');
            pg=permute(pg,izs);
            qax{1}=x;
            if any(m=='G') || czp>2 % image
                ngr=prod(szs(3:end));
                ncol=ceil(sqrt(ngr));
                nrow=ceil(ngr/ncol);
                npix=szs(1)*szs(2);
                ia=0;
                for i=1:ngr
                    subplot(nrow,ncol,i);
                    ib=ia+npix;
                    imagesc(qax{izs(1)},qax{izs(2)},reshape(pg(ia+1:ib),szs(1:2))');
                    axis 'xy'
                    if any(m=='c')
                        colorbar;
                    end
                    if nrow*ncol-i<ncol
                        xlabel(axlab(izs(1)));
                    end
                    if rem(i-1,ncol)==0
                        ylabel(axlab(izs(2)));
                    end
                    ia=ib;
                end
            else                    % graph
                plot(qax{izs(1)},reshape(permute(pg,izs),szs(1:2)),'-');
                xlabel(axlab{izs(1)});
            end
        else
            if any(m=='G')  % image
                imagesc(pg');
                axis 'xy'
                if any(m=='c')
                    colorbar;
                end
                xlabel('Input SNR (dB)');
                ylabel('Model Index');
            else            % graph
                if nx>=nq
                    plot(x,pg,'-');
                    xlabel('Input SNR (dB)');
                else
                    plot(1:nq,pg','-');
                    xlabel('Model Index');
                end
            end
        end
    end
end


