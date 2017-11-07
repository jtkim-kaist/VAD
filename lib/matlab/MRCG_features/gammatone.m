function r = gammatone(in, numChan, fRange, fs)
% Produce an array of filtered responses from a Gammatone filterbank.
% The first variable is required. 
% numChan: number of filter channels.
% fRange: frequency range.
% fs: sampling frequency.
% Written by ZZ Jin, adapted by DLW in Jan'07 and JF Woodruff in Nov'08

if nargin < 2
    numChan = 128;       % default number of filter channels in filterbank
end
if nargin < 3
    fRange = [80, 5000]; % default frequency range in Hz
end
if nargin < 4
    fs = 16000;     % default sampling frequency
end

filterOrder = 4;    % filter order
gL = 2048;          % gammatone filter length or 128 ms for 16 kHz sampling rate

sigLength = length(in);     % input signal length

phase(1:numChan) = zeros(numChan,1);        % initial phases
erb_b = hz2erb(fRange);       % upper and lower bound of ERB
erb = [erb_b(1):diff(erb_b)/(numChan-1):erb_b(2)];     % ERB segment
cf = erb2hz(erb);       % center frequency array indexed by channel
b = 1.019*24.7*(4.37*cf/1000+1);       % rate of decay or bandwidth

% Generating gammatone impulse responses with middle-ear gain normalization

gt = zeros(numChan,gL);  % Initialization
tmp_t = [1:gL]/fs;
for i = 1:numChan
    gain = 10^((loudness(cf(i))-60)/20)/3*(2*pi*b(i)/fs).^4;    % loudness-based gain adjustments
    gt(i,:) = gain*fs^3*tmp_t.^(filterOrder-1).*exp(-2*pi*b(i)*tmp_t).*cos(2*pi*cf(i)*tmp_t+phase(i));
end

sig = reshape(in,sigLength,1);      % convert input to column vector

% gammatone filtering using FFTFILT
r = fftfilt(gt',repmat(sig,1,numChan))';
