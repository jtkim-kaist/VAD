function ouotput = MRCG_features(sig, sampFreq)
% This function computes MRCG features

beta = 1000 / sqrt( sum(sig .^ 2) / length(sig) );
sig = sig .* beta;
sig = reshape(sig, length(sig), 1);
g = gammatone(sig, 64, [50 8000], sampFreq); % Gammatone filterbank responses

cochlea1 = log10(cochleagram(g,sampFreq*0.025,sampFreq*0.010));
cochlea2 = log10(cochleagram(g,sampFreq*0.200,sampFreq*0.010));
cochlea1 = cochlea1(:,:);
cochlea2 = cochlea2(:,:);

cochlea3  = get_avg(cochlea1,5,5);
cochlea4  = get_avg(cochlea1,11,11);
all_cochleas = [cochlea1; cochlea2; cochlea3; cochlea4];

del = deltas(all_cochleas);
ddel = deltas(deltas(all_cochleas,5),5);

ouotput = [all_cochleas;del;ddel];
