function loud=loudness(freq)
% Compute loudness level in Phons on the basis of equal-loudness functions.
% It accounts a middle ear effect and is used for frequency-dependent gain adjustments.
% This function uses linear interpolation of a lookup table to compute the loudness level, 
% in phons, of a pure tone of frequency freq using the reference curve for sound 
% pressure level dB. The equation is taken from section 4 of BS3383.
% 
% Written by ZZ Jin, and adapted by DLW in Jan'07

dB=60;

load('f_af_bf_cf.mat');
% Stores parameters of equal-loudness functions from BS3383,"Normal equal-loudness level
% contours for pure tones under free-field listening conditions", table 1.
% f (or ff) is the tone frequency, af and bf are frequency-dependent coefficients, and
% tf is the threshold sound pressure level of the tone, in dBs   

if (freq<20|freq>12500)
    return;
end
i=1;
while(ff(i)<freq)
    i=i+1;
end
afy=af(i-1)+(freq-ff(i-1))*(af(i)-af(i-1))/(ff(i)-ff(i-1));
bfy=bf(i-1)+(freq-ff(i-1))*(bf(i)-bf(i-1))/(ff(i)-ff(i-1));
cfy=cf(i-1)+(freq-ff(i-1))*(cf(i)-cf(i-1))/(ff(i)-ff(i-1));
loud=4.2+afy*(dB-cfy)/(1+bfy*(dB-cfy));