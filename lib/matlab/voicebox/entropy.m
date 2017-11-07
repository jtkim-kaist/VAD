function h=entropy(p,dim,cond,arg,step)
%ENTROPY calculates the entropy of discrete and sampled continuous distributions H=(P,DIM,COND,ARG,STEP)
%
%  Inputs:  P        is a vector or matrix of probabilities - one dimension per variable
%           DIM      lists dimensions along which to evaluate the entropy [default: 1st non singleton dimension]
%           COND     lists dimensions to use as conditional variables [default - none]
%           ARG      lists dimensions to use as parameters in the ouput [default - none]
%           STEP     for continuous distributions STEP gives the sample increment for each dimension of P
%                    if STEP is a scalar, the increment is assumed to be the same for each dimension
%
% Outputs:  H        is the entropy. It will have the same number of dimensions as the length of the ARG input.
%                    If the STEP argument is specified then this will be the differential entropy.
%
% Example: Suppose P(W,X,Y,Z) represents the joint probability of four correlated random variables
%
%               (a) H(W,X,Y,Z) = entropy(P,[1 2 3 4]). 
%               (b) H(W) = entropy(P), or equivalently entropy(P,1)
%               (c) H(W,Z | X,Y) = entropy(P,[1 4],[2 3])
%               (d) H(W | X, Z=z) = entropy(P,1,2,4); this is a function of z and will be a column vector
%
% As a special case, if the dimensions included in DIM are all singletons, the entries in P are treated
% as Bernoulli variable probabilities.

%	   Copyright (C) Mike Brookes 2006
%      Version: $Id: entropy.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<5
    stp=zeros(ndims(p),1);
else
    stp=repmat(step(1),ndims(p),1);
    stp(1:length(step))=step(:);
    stp=log2(stp);
end
if nargin<4
    arg=[];
else
    arg=arg(arg>0);
end
if nargin<3
    cond=[];
else
    cond=cond(cond>0);
end
if ~length(cond)
    s=size(p);
    if nargin<2
        dim=find(s>=min(2,max(s)));
        dim=dim(1);
    else
        dim=dim(dim>0);
    end
    st=prod(s);
    sd=prod(s(dim));
    sa=prod(s(arg));
    marg=1:length(s);
    marg(arg)=0;
    marg(dim)=0;
    marg=marg(marg>0);
    sm=st/sd/sa;
    if sm>1
        ip=[arg dim(:)' marg(:)'];
        sp=[s([arg dim(:)']) prod(s(marg))];
        q=sum(reshape(permute(p,[arg(:)' dim(:)' marg]),sa,sd,sm),3);
    else
        q=reshape(permute(p,[arg(:)' dim(:)' marg]),sa,sd);
    end
    if sd==1
        h=-log2(q+(q==0)).*q-log2(1-q+(q==1)).*(1-q);   % special treatment for bernoulli variables
    else
        sq=sum(q,2);
        h=sum(-log2(q+(q==0)).*q,2)./sq+log2(sq);
    end
    if length(arg)>1
        h=reshape(h,s(arg));
    end
    h=h+sum(stp(dim));
else
    % we could probably make this more efficient by avoiding the recursive call
    h=entropy(p,[dim(:); cond(:)],0,arg)-entropy(p,cond,0,arg);
end
