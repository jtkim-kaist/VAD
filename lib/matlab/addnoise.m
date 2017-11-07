function [ noisy_speech ] = addnoise( clean_speech, noise, SNR )

n = length(clean_speech);

% ensure the mean is zero
av = sum(noise)/n;
noise = noise - av;

% ensure the standard deviation is unity
var = sum(noise.^2)/n;
sd = sqrt(var);
noise = noise/sd;

% determine the power of the speech signal
var = sum(clean_speech.^2)/n;
ratio = 10^(SNR/10);
sd = sqrt(var/ratio);

% add noise to speech

noisy_speech = clean_speech + (noise*sd);

end

