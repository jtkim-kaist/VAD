function [y,ty]=correlogram(x,inc,nw,nlag,m,fs)
% make correlogram,
% Usage:
%        do_env=1; do_hp=1;                            % flags to control options
%        [b,a,fx,bx,gd]=gammabank(25,fs,'',[80 5000]); % determine filterbank
%        y=filterbank(b,a,s,gd);                       % filter signal s
%        if do_env
%            [bl,al]=butter(3,2*800/fs);
%            y=filter(bl,al,teager(y,1),[],1);           % low pass envelope @ 800 Hz
%        end
%        if do_hp
%            y=filter(fir1(round(16e-3*fs),2*64/fs,'high'),1,y,[],1);  % HP filter @ 64 Hz
%        end
%        correlogram(y,round(10e-3*fs),round(16e-3*fs),round(12e-3*fs),'',fs);
%
% Inputs:
%        x(*,chan)  is the output of a filterbank
%                   with one column per filter channel
%        inc        frame increment in samples
%        nw         window length in samples [or window function]
%        nlag       number of lags to calculate
%        m          mode string
%              [d = subtract DC component]
%              [n = unnormalized]
%              [v = variable analysis window proportional to lag]
%              [p = output the peaks only]
%               'h' = Hamming window
%        fs         sample freq (affects only plots)
%
% Outputs:
%        y(lag,chan,frame) is correlogram. Lags are 1:nlag samples
%        ty                is time of the window energy centre for each frame
%                            x(n) is at time n

memsize=voicebox('memsize');    % set memory size to use
[nx,nc]=size(x);  % number of sampes and channels
if nargin<6
    fs=1;
    if nargin<5
        m='h';
        if nargin<4
            nlag=[];
            if nargin<3
                nw=[];
                if nargin<2
                    inc=[];
                end
            end
        end
    end
end
if ~numel(inc)
    inc=128;
end
if ~numel(nw)
    nw=inc;
end
nwin=length(nw);
if nwin>1
    win=nw(:);
else
    nwin=nw;
    if any(m=='h')
        win=hamming(nwin);
    else
        win=ones(nwin,1);               % window function
    end
end
if ~numel(nlag)
    nlag=nwin;
end
nwl=nwin+nlag-1;
nt=pow2(nextpow2(nwl));  % transform length
nf=floor((nx-nwl+inc)/inc);  % number of frames
i1=repmat((1:nwl)',1,nc)+repmat(0:nx:nx*nc-1,nwl,1);
nb=min(nf,max(1,floor(memsize/(16*nwl*nc))));    % chunk size for calculating
nl=ceil(nf/nb);                  % number of chunks
jx0=nf-(nl-1)*nb;                % size of first chunk in frames
wincg=(1:nwin)*win.^2/(win'*win);
fwin=conj(fft(win,nt,1)); % fft of window
y=zeros(nlag,nc,nf);
% first do partial chunk
jx=jx0;
x2=zeros(nwl,nc*jx);
x2(:)=x(repmat(i1(:),1,jx)+repmat((0:jx-1)*inc,nwl*nc,1));
v=ifft(conj(fft(x2(1:nwin,:),nt,1)).*fft(x2,nt,1));
w=real(ifft(fwin(:,ones(1,nc*jx)).*fft(x2.*conj(x2),nt,1)));
w=sqrt(w(1:nlag,:).*w(ones(nlag,1),:));
if isreal(x)
    y(:,:,1:jx)=reshape(real(v(1:nlag,:))./w,nlag,nc,jx);
else
    y(:,:,1:jx)=reshape(conj(v(1:nlag,:))./w,nlag,nc,jx);
end
% now do remaining chunks
x2=zeros(nwl,nc*nb);
for il=2:nl
    ix=jx+1;            % first frame in this chunk
    jx=jx+nb;           % last frame in this chunk
    x2(:)=x(repmat(i1(:),1,nb)+repmat((ix-1:jx-1)*inc,nwl*nc,1));
    v=ifft(conj(fft(x2(1:nwin,:),nt,1)).*fft(x2,nt,1));
    w=real(ifft(fwin(:,ones(1,nc*nb)).*fft(x2.*conj(x2),nt,1)));
    w=sqrt(w(1:nlag,:).*w(ones(nlag,1),:));
    if isreal(x)
        y(:,:,ix:jx)=reshape(real(v(1:nlag,:))./w,nlag,nc,nb);
    else
        y(:,:,ix:jx)=reshape(conj(v(1:nlag,:))./w,nlag,nc,nb);
    end
end
ty=(0:nf-1)'*inc+wincg;       % calculate times of window centres
if ~nargout
    imagesc(ty/fs,(1:nlag)/fs,squeeze(mean(y,2)));
    if nargin<6
        us='samp';
    else
        us='s';
    end
    xlabel(['Time (' xticksi us ')']);
    ylabel(['Lag (' yticksi us ')']);
    axis 'xy';
    title('Summary Correlogram');
end


