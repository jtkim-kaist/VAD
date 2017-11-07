function [xx,ii,m,v]=psycestu(iq,x,r,xp)
% psycestu estimate unimodal psychometric function
%
% Usage: [xx,ii,m,v]=psycestu(-n,p,q,xp) % initialize n models
%        [xx,ii,m,v]=psycestu(i,x,r)     % supply a trial result to psycest
%                    psycestu(i)         % plot pdf of model i
%              [p,q]=psychestu(0)        % output model parameters (or print them if no outputs)
%
% Inputs:
%         -n        minus the number of models
%          p(:,n)   parameters for each model
%                      1  miss [0.04]
%                      2  guess [0.1]
%                      3,4  SNR at peak: [min max] [-20 20]
%                      5,6  normalized semi-width: [min max] [0 20]
%                      7,8  low side slope: [min max] [0.03 0.3]
%                      9,10  high side slope: [min max] [0.03 0.3]
%          q(:)     parameters common to all models (vector or struct)
%                      1  q.nxslh  size of pdf array: [nx ns nl nh] [20 21 19 18]
%                      2  q.nh  number of probe SNR values to evaluate [30]
%                      3  q.cs  weighting of pdf factors [1 1 1 1]
%                      5  q.dh  minimum step size in dB for probe SNRs [0.2]
%                      6  q.sl  min slopes threshold [0.02]
%                      7  q.kp  number of std deviations of the pdfs to keep [4]
%                      8  q.hg  amount to grow expected gains in ni trials [1.3]
%          xp{n}(:) list of available probe SNR values
%          i        model probed
%          x        probe SNR value used
%          r        response: 0=wrong, 1=right.
%
% Outputs:
%          xx       recommended probe SNR
%          ii       recommended model to probe next
%          m(4,n,3) estimated srt and slope of all models
%                   m(:,:,1:3) are respectively the mean, mode (MAP) and marginal mode estimates
%          v(4,n)   estimated diagonal covariance matrix entries:
persistent pp qq fl fh fxs tx ts tl th
if iq<0
    if iq~=-1
        error('Cannot yet have multiple models');
    end
    pp=x;
    qq=r;
    pp(7:2:9)=max(pp(7:2:9),qq.sl);  % enforce the minimum slope
end
nxslh=qq.nxslh;  % make local copies of some useful values
nx=nxslh(1);
ns=nxslh(2);
nl=nxslh(3);
nh=nxslh(4);
la=pp(1);
gu=pp(2);
if iq<0
    % reserve space for pdfs
    fl=ones(nx,nl,ns);
    fh=ones(nx,nh,ns);
    fxs=ones(nx,1,ns); % marginal x-s distribution
    % define initial ranges
    tx=linspace(pp(3),pp(4),nx)';
    ts=linspace(pp(5),pp(6),ns)';
    tl=linspace(log(pp(7)),log(pp(8)),nl)';
    th=linspace(log(pp(9)),log(pp(10)),nh)';
elseif iq>0 && nargin==3
    if iq~=1
        error('Cannot yet have multiple models');
    end
    r0=r==0;
    % Update the pdf arrays
    mskx=tx>x; % find which tx values are >x
    nxm=sum(mskx);
    ets=reshape(exp(-ts),[1 1 ns]);
    sd0=1+ets.^2;
    sm0=2*ets;
    if any(mskx)
        fl(mskx,:,:)=fl(mskx,:,:).*(r0+(1-2*r0)*(gu+(1-gu-la)*(repmat(sd0,[nxm,nl,1])+repmat(sm0,[nxm,nl,1]).*repmat(cosh((tx(mskx)-x)*exp(tl')),[1 1 ns])).^(-1)));
    end
    if ~all(mskx)
        fh(~mskx,:,:)=fh(~mskx,:,:).*(r0+(1-2*r0)*(gu+(1-gu-la)*(repmat(sd0,[nx-nxm,nh,1])+repmat(sm0,[nx-nxm,nh,1]).*repmat(cosh((x-tx(~mskx))*exp(th')),[1 1 ns])).^(-1)));
    end
else % iq=0
    xx=pp;
    ii=qq;
end

% Normalize the pdf arrays

sfl=sum(fl,2);  % sum over slope variable
fl=fl./sfl(:,ones(nl,1),:); % normalize each x-s plane
sfh=sum(fh,2);
fh=fh./sfh(:,ones(nh,1),:); % normalize each x-s plane
fxs=fxs.*sfl.*sfh;          % Marginal x-s distribution
fxs=fxs/sum(fxs(:));        % Normalize to 1

% calculate marginal distributions

px=squeeze(sum(fxs,3));
ps=squeeze(sum(fxs,1));
pl=reshape(permute(fl,[2 1 3]),nl,nx*ns)*fxs(:);
ph=reshape(permute(fh,[2 1 3]),nh,nx*ns)*fxs(:);

% calculate means, modes and entropies

m=zeros(4,1,3);
m(:,1,1)=[tx'*px ts'*ps tl'*pl th'*ph]';
vraw=[tx'.^2*px ts'.^2*ps tl'.^2*pl th'.^2*ph]'-m(:,1,1).^2;
% h=[(tx(1)-tx(2))*px'*log(px) (ts(1)-ts(2))*ps'*log(ps) (tl(1)-tl(2))*pl'*log(pl) (th(1)-th(2))*ph'*log(ph)]';
% calculate joint mode
[flm,iflm]=max(fl,[],2);
[fhm,ifhm]=max(fh,[],2);
fxsm=fxs.*flm.*fhm;
[pxpk,i]=max(fxsm(:)); % find global mode
j=1+floor((i-1)/nx);
i=i-nx*(j-1);
jl=iflm(i,1,j);
jh=ifhm(i,1,j);
% could use quadratic interpolation but it seems a bit flaky for 4-D
% if all([i j jl jh]>1) && all([i j jl jh]<[nx ns nl nh])
%     [dum,ij]=quadpeak(repmat(fxs(i-1:i+1,1,j-1:j+1),[1 3 1 3]).*repmat(fl(i-1:i+1,jl-1:jl+1,j-1:j+1),[1 1 1 3]).*permute(repmat(fh(i-1:i+1,jh-1:jh+1,j-1:j+1),[1 1 1 3]),[1 4 3 2]));
% end
m(:,1,2)=[tx(i) ts(j) tl(jl) th(jh)]'; % replicate mean for now
% calculate marginal modes
[pxpk,i]=max(px);                          % marginal mode in x
if i>1 && i<nx                           % use quadratic interpolation if possible
    [dum,j]=quadpeak(px(i-1:i+1));
    i=i+j-2;
end
xmm=(2-i)*tx(1)+(i-1)*tx(2);             % marginal mode in x

[pxpk,i]=max(ps);                          % marginal mode in sw
if i>1 && i<ns                           % use quadratic interpolation if possible
    [dum,j]=quadpeak(ps(i-1:i+1));
    i=i+j-2;
end
smm=(2-i)*ts(1)+(i-1)*ts(2);             % marginal mode in sw

[pxpk,i]=max(pl);                          % marginal mode in l
if i>1 && i<nl                           % use quadratic interpolation if possible
    [dum,j]=quadpeak(pl(i-1:i+1));
    i=i+j-2;
end
lmm=(2-i)*tl(1)+(i-1)*tl(2);             % marginal mode in l

[pxpk,i]=max(ph);                          % marginal mode in h
if i>1 && i<nh                           % use quadratic interpolation if possible
    [dum,j]=quadpeak(ph(i-1:i+1));
    i=i+j-2;
end
hmm=(2-i)*th(1)+(i-1)*th(2);             % marginal mode in h
m(:,1,3)=[xmm smm lmm hmm]';
m(3:4,:,:)=exp(m(3:4,:,:));              % convert log slopes to slopes
mfact=(m(3,1,:)+m(4,1,:))./(m(3,1,:).*m(4,1,:));
m(2,1,:)=m(2,1,:).*mfact;  % convert normalized semi-width to width
v=[vraw(1); vraw(2)*mfact(1).^2; vraw(3:4).*m(3:4,1,1).^2]; % calculate variances

% figure(21); plot(tx,px,m([1 1]),[0 1.05*max(px)]); title('SNR at peak');
% figure(22); plot(ts,ps,m([2 2]),[0 1.05*max(ps)]); title('semi-width');
% figure(23); plot(tl,pl,log(m([3 3])),[0 1.05*max(pl)]); title('lower log slope');
% figure(24); plot(th,ph,log(m([4 4])),[0 1.05*max(ph)]); title('upper log slope');
if ~nargout && iq==1  % plot pdf
    subplot(121)
    imagesc(tx,ts,squeeze(fxs)');
    axis 'xy'
    xlabel('Peak position (dB)');
    ylabel('Normalized semi-width (dB)');
    subplot(122)
    imagesc(th,tl,reshape(permute(fl,[2 1 3]),nl,nx*ns)*reshape(permute(fh.*fxs(:,ones(nh,1),:),[1 3 2]),nx*ns,nh));
        axis 'xy'
        xlabel('Ln down slope (ln prob/dB)');
         ylabel('Ln up slope (ln prob/dB)');
end
% now determine the next recommended probe SNR

if iq~=0
xt=tx;  % for now we just try all the tx values
nxt=numel(xt);
ets=reshape(exp(-ts),[1 1 ns]);
sd0=1+ets.^2;
sm0=2*ets;
hh=zeros(nxt,1);  % store for entropies
for i=1:nxt
    y=xt(i);        % y = probe value
    mskx=tx>y; % find which tx values are >y
    nxm=sum(mskx);
    flm=gu+(1-gu-la)*(repmat(sd0,[nxm,nl,1])+repmat(sm0,[nxm,nl,1]).*repmat(cosh((tx(mskx)-y)*exp(tl')),[1 1 ns])).^(-1);
    fhm=gu+(1-gu-la)*(repmat(sd0,[nx-nxm,nh,1])+repmat(sm0,[nx-nxm,nh,1]).*repmat(cosh((y-tx(~mskx))*exp(th')),[1 1 ns])).^(-1);
    % calculate probs conditional on r=1
    fl1=fl;
    fh1=fh;
    fl1(mskx,:,:)=fl1(mskx,:,:).*flm;
    fh1(~mskx,:,:)=fh1(~mskx,:,:).*fhm;
    sfl1=sum(fl1,2);  % sum over slope variable
    fl1=fl1./sfl1(:,ones(nl,1),:); % normalize each x-s plane
    sfh1=sum(fh1,2);
    fh1=fh1./sfh1(:,ones(nh,1),:); % normalize each x-s plane
    fxs1=fxs.*sfl1.*sfh1;          % Marginal x-s distribution
    pr=sum(fxs1(:)); % P(r=1 | y)
    fxs1=fxs1/pr;        % Normalize to 1
    px1=squeeze(sum(fxs1,3));
    ps1=squeeze(sum(fxs1,1));
    pl1=reshape(permute(fl1,[2 1 3]),nl,nx*ns)*fxs1(:);
    ph1=reshape(permute(fh1,[2 1 3]),nh,nx*ns)*fxs1(:);
    % calculate probs conditional on r=0
    fl0=fl;
    fh0=fh;
    fl0(mskx,:,:)=fl0(mskx,:,:).*(1-flm);
    fh0(~mskx,:,:)=fh0(~mskx,:,:).*(1-fhm);
    sfl0=sum(fl0,2);  % sum over slope variable
    fl0=fl0./sfl0(:,ones(nl,1),:); % normalize each x-s plane
    sfh0=sum(fh0,2);
    fh0=fh0./sfh0(:,ones(nh,1),:); % normalize each x-s plane
    fxs0=fxs.*sfl0.*sfh0;          % Marginal x-s distribution
    fxs0=fxs0/(1-pr);        % Normalize to 1
    px0=squeeze(sum(fxs0,3));
    ps0=squeeze(sum(fxs0,1));
    pl0=reshape(permute(fl0,[2 1 3]),nl,nx*ns)*fxs0(:);
    ph0=reshape(permute(fh0,[2 1 3]),nh,nx*ns)*fxs0(:);
    % now calculate the entropies
    h1=[(tx(1)-tx(2))*px1'*log(px1) (ts(1)-ts(2))*ps1'*log(ps1) (tl(1)-tl(2))*pl1'*log(pl1) (th(1)-th(2))*ph1'*log(ph1)]';
    h0=[(tx(1)-tx(2))*px0'*log(px0) (ts(1)-ts(2))*ps0'*log(ps0) (tl(1)-tl(2))*pl0'*log(pl0) (th(1)-th(2))*ph0'*log(ph0)]';
    hh(i)=qq.cs*(pr*h1+(1-pr)*h0);
end
[hmin,ih]=min(hh);
xx=xt(ih);
ii=1;
end

% now rescale the pdf arrays

% sraw=sqrt(vraw);  % std deviations
% clim=[tx(1) tx(end); ts(1) ts(end); tl(1) tl(end); th(1) th(end)]  % current axis limits
% minlim=clim*[3 1; 0 2]/3; % always keep at least the middle third of the existing array
% maxlim=clim*[2 0; 1 3]/3;
% pra=min(max(repmat(m(:,1,1),1,2)+qq.kp*sraw*[-1 1],minlim),maxlim);

plow=max(min(0.1*gu,0.001),1e-6);

pcx=cumsum(px);
tmp=linspace(tx(max(find(pcx>plow,1)-1,1)),tx(find(pcx>(1-plow),1)),nx)';
fxs=linres(fxs,1,tx,tmp);
fl=linres(fl,1,tx,tmp);
fh=linres(fh,1,tx,tmp);
tx=tmp;
%
pcx=cumsum(ps);
tmp=linspace(ts(max(find(pcx>plow,1)-1,1)),ts(find(pcx>(1-plow),1)),ns)';
fxs=linres(fxs,3,ts,tmp);
fl=linres(fl,3,ts,tmp);
fh=linres(fh,3,ts,tmp);
ts=tmp;
%
% For now, we do not update the slope distributions
%
% tmp=linspace(pra(3,1),pra(3,2),nl)';
% fl=linres(fl,2,tl,tmp);
% tl=tmp;
%
% tmp=linspace(pra(4,1),pra(4,2),nh)';
% fh=linres(fh,2,th,tmp);
% th=tmp;

function y=linres(x,d,a,b)
% linear resample x along dimension d from grid a to b
sz=size(x);
n=sz(d);
p=[d:numel(sz) 1:d-1];
sz2=prod(sz)/n;
z=reshape(permute(x,p),n,sz2);  % put the target dimension in column 1
na=numel(a);
nb=numel(b);
[xx,ix]=sort([a(:); b(:)]);
jx=zeros(1,na+nb);
jx(ix)=(1:na+nb);
% ja=jx(1:na)-(1:na); % last new point < each old point
jb=jx(na+1:end)-(1:nb); % last old point < each new point
y=zeros(nb,size(z,2));
y(jb==0,:)=z(ones(sum(jb==0),1),:);  % replicated entries
y(jb==na,:)=z(na*ones(sum(jb==na),1),:);
mk=(jb>0) & (jb<na);  % interpolation mask
jmk=jb(mk);
f=((b(mk)-a(jmk))./(a(1+jmk)-a(jmk)));
fc=1-f;
y(mk,:)=z(jmk,:).*fc(:,ones(1,sz2))+z(1+jmk,:).*f(:,ones(1,sz2));
sz(d)=nb;
y=ipermute(reshape(y,sz(p)),p);






