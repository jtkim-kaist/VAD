function [features, gfilters] = FE_GBF(signal, fs , nb_mod_freq, mvnon, gtfeat)
%usage: [features, gfilters] = gbfb_fe(signal, fs)
%   signal:   waveform signal
%   fs:       sampling rate in Hz
% 
%   features: Gabor filter bank (GBFB) features 
%   gfilters: Gabor filter bank filters
%
% - Gabor Filter Bank Feature Extraction v2.0 -
%
% Autor    : Marc René Schädler
% Email    : marc.r.schaedler@uni-oldenburg.de
% Institute: Medical Physics / Carl-von-Ossietzky University Oldenburg, Germany
% 
%-----------------------------------------------------------------------------
%
% Licensing of the feature extraction code: Dual-License
% 
% The Gabor Filterbank (GBFB) feature extraction code is licensed under both 
% General Public License (GPL) version 3 and a proprietary license that can be 
% arranged with us. In practical sense, this means:
% 
% - If you are developing Open Source Software (OSS) based on the GBFB code, 
%   chances are you will be able to use it freely under GPL. But please double check 
%   http://www.gnu.org/licenses/license-list.html for OSS license compatibility with GPL
% - Alternatively, if you are unable to release your application as Open Source Software, 
%   you may arrange alternative licensing with us. Just send your inquiry to the 
%   author to discuss this option.
% 
%-----------------------------------------------------------------------------
% 
% Release Notes:
% This script is written to produce the results published in [1]
% For an overview of current publications / further developments of this 
% frontend, visit http://medi.uni-oldenburg.de/GBFB
% 
% [1] M.R. Schädler, B.T. Meyer, B. Kollmeier
% "Spectro-temporal modulation subspace-spanning filter bank features 
% for robust automatic speech recognition ", J. Acoust. Soc. Am. Volume 131, 
% Issue 5, pp. 4134-4151 (2012)
%
% Paper URL: http://link.aip.org/link/?JAS/131/4134
% Paper DOI: 10.1121/1.3699200
%
%-----------------------------------------------------------------------------
% Adapted by:  Maarten Van Segbroeck
% Institute: University of Southern California, SAIL, U.S
% Mail: maarten@sipi.usc.edu
% 

%% Default settings

% Filter bank settings [spectral temporal]
omega_max   = [pi/2 pi/2]; 	% radian
size_max    = [3*23   40]; 	% bands, frames
nu          = [ 3.5  3.5]; 	% half-waves under envelope
distance    = [ 0.3  0.2]; 	% controls the spacing of filters

%% ADAPTION [MVS]
if ~exist('nb_mod_freq', 'var')
    numGFCCs = 2;  % selection of first number of temporal modulation frequencies
end
if ~exist('mvnon', 'var')
    mvnon = true;  % selection of first number of temporal modulation frequencies
end


%% Set up Gabor filter bank with given parameters

% Calculate center modulation frequencies.
[omega_n, omega_k] = gfbank_calcaxis(omega_max, size_max, nu, distance);
%% ADAPTION [MVS] - SELECT ONLY FILTER WITH LOW TEMPORAL MODULATION
omega_n=omega_n(1:nb_mod_freq);

omega_n_num = length(omega_n);
omega_k_num = length(omega_k);

% Generate filters for all pairs of spectral and temporal modulation
% frequencies except for the redundant ones.
gfilters = cell(omega_k_num,omega_n_num);
for i=1:omega_k_num
    for j=1:omega_n_num
        if ~(omega_k(i) < 0 && omega_n(j) == 0)
            gfilters{i,j} = gfilter_gen(omega_k(i), omega_n(j), nu(1), nu(2), size_max(1), size_max(2));
        end
    end
end

%% ADAPTION [MVS] - PAD WITH ZEROS
paddingw = floor(min(0.2*(fs/8000)*fs,floor((length(signal)-1)*100/fs)*fs/100)); % 200 ms padding window
paddingf = paddingw/fs*100; % frameshift: 10ms
signal=reshape(signal,1,length(signal)); % make data consistent
signal = [signal(1:paddingw) signal signal(end-paddingw:end)];

%% ADAPTION [MVS] -  USE ALTERNATIVE MEL-SPECTROGRAM
if ~exist('gtfeat', 'var')
    mel_spec_log = FE_MEL(signal,fs);
else
    mel_spec_log = gtfeat;
end

%% Filter mel spectrogram with filter bank filters and select representative channels.
gfilters_output = cell(omega_k_num,omega_n_num);
for i=1:omega_k_num
    for j=1:omega_n_num
        gfilter = gfilters{i,j};
        if ~isempty(gfilter)
            % Filter mel spectrogram with Gabor filter.
            mel_spec_log_filtered = gfilter_filter(gfilter, mel_spec_log);
            % Select representative channels from filtered Mel-spectrogram.
            gfilters_output{i,j} = gfilter_rep(gfilter, mel_spec_log_filtered);
        end
    end
end
features = cell2mat(reshape(gfilters_output,[],1));

% Use the real part of the filter output
features = real(features);

%% ADAPTION [MVS] - REMOVE PADDING
if ~exist('gtfeat', 'var')
    features = features(:,paddingf+1:end-paddingf);
end

%% ADAPTION [MVS] - APPLY MEAN AND VARIANCE NORMALIZATION (DEFAULT ON)
if mvnon
features = mvn(features);
end

end


function [omega_n, omega_k] = gfbank_calcaxis(omega_max, size_max, nu, distance)
% Calculates the modulation center frequencies iteratively.
% Termination condition for iteration is reaching omega_min, which is
% derived from size_max.
omega_min = (pi .* nu) ./ size_max;

% Eq. (2b)
c = distance .* 8 ./ nu;
% Second factor of Eq. (2a)
space_n = (1 + c(2)./2) ./ (1 - c(2)./2);
count_n = 0;
omega_n(1) = omega_max(2);
% Iterate starting with omega_max in spectral dimension
while omega_n/space_n > omega_min(2)
    omega_n(1+count_n) = omega_max(2)/space_n.^count_n;
    count_n = count_n + 1;
end
omega_n = fliplr(omega_n);
% Add DC
omega_n = [0,omega_n];
% Second factor of Eq. (2a)
space_k = (1 + c(1)./2) ./ (1 - c(1)./2);
count_k = 0;
omega_k(1) = omega_max(1);
% Iterate starting with omega_max in temporal dimension
while omega_k/space_k > omega_min(1)
    omega_k(1+count_k) = omega_max(1)/space_k.^count_k;
    count_k = count_k + 1;
end
% Add DC and negative MFs for spectro-temporal opposite 
% filters (upward/downward)
omega_k = [-omega_k,0,fliplr(omega_k)];
end


function gfilter = gfilter_gen(omega_k, omega_n, nu_k, nu_n, size_max_k, size_max_n)
% Generates a gabor filter function with:
%  omega_k       spectral mod. freq. in rad
%  omega_n       temporal mod. freq. in rad
%  nu_k          number of half waves unter the envelope in spectral dim.
%  nu_n          number of half waves unter the envelope in temporal dim.
%  size_max_k    max. allowed extension in spectral dimension
%  size_max_n    max. allowed extension in temporal dimension

% Calculate windows width.
w_n = 2*pi / abs(omega_n) * nu_n / 2;
w_k = 2*pi / abs(omega_k) * nu_k / 2;

% If the size exceeds the max. allowed extension in a dimension set the
% corresponding mod. freq. to zero.
if w_n > size_max_n
    w_n = size_max_n;
    omega_n = 0;
end
if w_k > size_max_k
    w_k = size_max_k;
    omega_k = 0;
end

% Separable hanning envelope, cf. Eq. (1c).
env_n = hann_win(w_n-1);
env_k = hann_win(w_k-1);
envelope = env_k * env_n.';
[win_size_k, win_size_n] = size(envelope);

% Sinusoid carrier, cf. Eq. (1c).
n_0 = (win_size_n+1) / 2;
k_0 = (win_size_k+1) / 2;
[n,k] = meshgrid (1:win_size_n, 1:win_size_k);
sinusoid = exp(1i*omega_n*(n - n_0) + 1i*omega_k*(k - k_0));

% Eq. 1c
gfilter  = envelope .* sinusoid;

% Compensate the DC part by subtracting an appropiate part
% of the envelope if filter is not the DC filter.
envelope_mean = mean(mean(envelope));
gfilter_mean = mean(mean(gfilter));
if (omega_n ~= 0) || (omega_k ~= 0)
    gfilter = gfilter - envelope./envelope_mean .* gfilter_mean;
else
    % Add an imaginary part to DC filter for a fair real/imag comparison.
    gfilter = gfilter + 1i*gfilter;
end
% Normalize filter to have gains <= 1.
gfilter = gfilter ./ max(max(abs(fft2(gfilter))));

end


function log_mel_spec_filt = gfilter_filter(gfilter, log_mel_spec)
% Applies the filtering with a 2D Gabor filter to log_mel_spec
% This includes the special treatment of filters that do not lie fully
% inside the spectrogram
if any(any(gfilter < 0))
    % Compare this code to the compensation for the DC part in the
    % 'gfilter_gen' function. This is an online version of it removing the
    % DC part of the filters by subtracting an appropriate part of the
    % filters' envelope.
    gfilter_abs_norm = abs(gfilter) ./ sum(sum(abs(gfilter)));
    gfilter_dc_map = fftconv2(ones(size(log_mel_spec)), gfilter,'same');
    env_dc_map = fftconv2(ones(size(log_mel_spec)), gfilter_abs_norm,'same');
    dc_map = fftconv2(log_mel_spec, gfilter_abs_norm,'same') ./ env_dc_map .* gfilter_dc_map;
else
    % Dont' remove the DC part if it is the DC filter.
    dc_map = 0;
end
% Filter log_mel_spec with the 2d Gabor filter and remove the DC parts.
log_mel_spec_filt = fftconv2(log_mel_spec, gfilter,'same') - dc_map;
end


function mel_spec_rep = gfilter_rep(gfilter, mel_spec)
% Selects the center channel by choosing k_offset and those with k_factor
% channels distance to it in spectral dimension where k_factor is approx.
% 1/4 of the filters extension in the spectral dimension.
k_factor = floor(1/4 * size(gfilter,1));
if k_factor < 1
    k_factor = 1;
end
k_offset = mod(floor(size(mel_spec,1)/2),k_factor);
k_idx = (1+k_offset):k_factor:size(mel_spec,1);
% Apply subsampling.
mel_spec_rep = mel_spec(k_idx,:);
end

function featmvn = mvn(feat)
% Mean and Variance Normalization
feat_mean = mean(feat,2);
feat_std = sqrt(var(feat,1,2));
featmvn = (feat - repmat(feat_mean,1,size(feat,2))) ./ repmat(feat_std,1,size(feat,2));
end

function window_function = hann_win(width)
% A hanning window function that accepts non-integer width and always
% returns a symmetric window with an odd number of samples.
x_center = 0.5;
x_values = [fliplr((x_center-1/(width+1)):-1/(width+1):0), x_center:1/(width+1):1]';
valid_values_mask = (x_values > 0) & (x_values < 1);
window_function =  0.5*(1-(cos(2*pi*x_values(valid_values_mask))));
end

function out = fftconv2(in1, in2, shape)
% 2D convolution in terms of the 2D FFT that substitutes conv2(in1, in2, shape).
size_y = size(in1,1)+size(in2,1)-1;
size_x = size(in1,2)+size(in2,2)-1;
fft_size_x = 2.^ceil(log2(size_x));
fft_size_y = 2.^ceil(log2(size_y));
in1_fft = fft2(in1,fft_size_y,fft_size_x);
in2_fft = fft2(in2,fft_size_y,fft_size_x);
out_fft = in1_fft .* in2_fft;
out_padd = ifft2(out_fft);
out_padd = out_padd(1:size_y,1:size_x);
switch shape
    case 'same'
        y_offset = floor(size(in2,1)/2);
        x_offset = floor(size(in2,2)/2);
        out = out_padd(1+y_offset:size(in1,1)+y_offset,1+x_offset:size(in1,2)+x_offset);
    case 'full'
        out = out_padd;
end
end

% ------------------------------------------------------------------------
% Copyright (C) 2012 University of Southern California, SAIL, U.S.
% Author: Maarten Van Segbroeck
% Mail: maarten@sipi.usc.edu
% ------------------------------------------------------------------------

%% LOG-MEL SCALED SPECTRAL FEATURES
function [lp,logE]=FE_MEL_RATS(sam,fs,tMelMat)
% --IN--
% sam: signal
% fs: sampling frequency
% tMelMat: MEL matrix (optional)
% --OUT--
% lp: log-mel-power of sam
% logE: log-energy of sam

if ~exist('fs', 'var')
    fs = 16000;
end

FrameLen=0.02*fs;
FrameShift=0.01*fs;
Efloor=exp(-50);
MELfloor=exp(-50);
PreEmp=0.97;
Nfft=2^ceil(log(FrameLen)/log(2));
NbCh=floor(hz2meldm(fs/2))-1; % number of MEL filter banks 
if ~exist('tMelMat', 'var')
    tMelMat=mel_matssp(fs/2,NbCh,Nfft)';
end

% truncate
NbFr=floor( (length(sam)-FrameLen+FrameShift)/FrameShift);
sam=sam(1:NbFr*FrameShift+FrameLen-FrameShift);

% DC removal
sam=filter([1 -1],[1 -0.999],[0 reshape(sam,1,length(sam))]);

% framing
ind1=(1:FrameShift:NbFr*FrameShift)-1;
ind2=(1:FrameLen)';
x=sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr)));

% logE
logE=log(max(sum(x.^2,1),Efloor));

% preemphasis & windowing
T=length(sam);
sam=[0 sam(2:T)-PreEmp*sam(1:T-1)];
win=hamming(FrameLen);
xNoWin=sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr)));
xWin=win(:,ones(1,NbFr)).*xNoWin;

% logMel
Nfft2=2*Nfft; % focus on 0-4000Hz band
X=fft(xWin,Nfft2);
X(1,:)=0;
Xmag=abs(X(1:Nfft/2+1,:));
lp=log(max(tMelMat*Xmag,MELfloor));

end

%% MEL-SCALED FILTERBANK
function [M,fCen]=mel_matssp(fs,NbCh,Nfft)
% --IN--
% fs: sampling frequency
% NbCh: number of MEL filters
% Nfft: number of FFT points 
% --OUT--
% M: Nfft/2+1-by-NbCh matrix of triangular weights
% fcen: center frequencies of triangular weights

fCen=meldm2hz(1:NbCh);
fLow=meldm2hz(0:NbCh-1);
fHigh=meldm2hz(2:NbCh+1);

flsamp = max(hz2spc(fLow,fs,Nfft/2+1),0);
i=find(spc2hz(flsamp,fs,Nfft/2+1)<fLow);
flsamp(i)=flsamp(i)+1;

fhsamp = min(hz2spc(fHigh,fs,Nfft/2+1),Nfft/2);
i=find(spc2hz(fhsamp,fs,Nfft/2+1)>fHigh);
fhsamp(i)=fhsamp(i)-1;

TotLen=fhsamp-fhsamp+1;
M=sparse([],[],[],Nfft/2+1,NbCh,sum(TotLen));
for k=1:NbCh,
  f=spc2hz(flsamp(k):fhsamp(k),fs,Nfft/2+1);
  w=zeros(size(f));
  sel=f<=fCen(k);
  w(sel)=max(0,1+(f(sel)-fCen(k))/(fCen(k)-fLow(k)));
  sel=f>fCen(k);
  w(sel)=max(0,(f(sel)-fHigh(k))/(fCen(k)-fHigh(k)));
  M(flsamp(k)+1:fhsamp(k)+1,k)=w';
end

end
