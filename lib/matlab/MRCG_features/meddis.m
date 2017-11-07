function hc = meddis(r, fs)
% Produce auditory nerve response from output of a Gammatone filterbank.
% This function uses the Meddis hair cell model that transduces from filter
% response to auditory nerve firing activity (the same as inner hair cell
% response).
% The first augument is required.
% fs: sampling frequency.
% Written by ZZ Jin, and adapted by DLW in Jan'07

if nargin < 2
    fs = 16000;     % default sampling frequency
end

[numChan,sigLength] = size(r);     % number of channels and input signal length

% hair cell parameters
med_y=5.05;
med_g=2000;
med_l=2500;
med_r=6580;
med_x=66.31;
med_a=3.0;
med_b=200;
med_h=48000;
med_m=1;

% initialize inner hair cells
ymdt=med_y*med_m/fs;
xdt=med_x/fs;
ydt=med_y/fs;
lplusrdt=(med_l+med_r)/fs;
rdt=med_r/fs;
gdt=med_g/fs;
hdt=med_h;

% inner hair cell transduction
hc=zeros(numChan,sigLength);
for i=1:numChan
    kt=med_g*med_a/(med_a+med_b);
    hair_c=med_m*med_y*kt/(med_l*kt+med_y*(med_l+med_r));
    hair_q=hair_c*(med_l+med_r)/kt;
    hair_w=hair_c*med_r/med_x;

    for j=1:sigLength
        if (r(i,j)+med_a)>0
            kt=gdt*(r(i,j)+med_a)/(r(i,j)+med_a+med_b);
        else
            kt=0;
        end
        if hair_q<med_m
            replenish=ymdt-ydt*hair_q;
        else
            replenish=0;
        end
        eject=kt*hair_q;
        reuptakeandloss=lplusrdt*hair_c;
        reuptake=rdt*hair_c;
        reprocess=xdt*hair_w;
        hair_q=max(hair_q+replenish-eject+reprocess,0);
        hair_c=max(hair_c+eject-reuptakeandloss,0);
        hair_w=max(hair_w+reuptake-reprocess,0);
        hc(i,j)=hair_c*hdt;
    end
end