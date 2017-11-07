function [wts,cfreqs] = gammatone_matrix(nfft, sr, nfilts, width, minfreq, maxfreq, maxlen)
% [wts,cfreqa] = fft2gammatonemx(nfft, sr, nfilts, width, minfreq, maxfreq, maxlen)
%      Generate a matrix of weights to combine FFT bins into
%      Gammatone bins.  nfft defines the source FFT size at
%      sampling rate sr.  Optional nfilts specifies the number of
%      output bands required (default 64), and width is the
%      constant width of each band in Bark (default 1).
%      minfreq, maxfreq specify range covered in Hz (100, sr/2).
%      While wts has nfft columns, the second half are all zero. 
%      Hence, aud spectrum is
%      fft2gammatonemx(nfft,sr)*abs(fft(xincols,nfft));
%      maxlen truncates the rows to this many bins.
%      cfreqs returns the actual center frequencies of each
%      gammatone band in Hz.
%
% 2004-09-05  Dan Ellis dpwe@ee.columbia.edu  based on rastamat/audspec.m
% Last updated: $Date: 2009/02/22 02:29:25 $

if nargin < 2;    sr = 16000; end
if nargin < 3;    nfilts = 64; end
if nargin < 4;    width = 0.5; end
if nargin < 5;    minfreq = 50; end
if nargin < 6;    maxfreq = sr/2; end
if nargin < 7;    maxlen = nfft/2; end

wts = zeros(nfilts, nfft);

% after Slaney's MakeERBFilters
EarQ = 9.26449;
minBW = 24.7;
order = 1;

cfreqs = -(EarQ*minBW) + exp((1:nfilts)'*(-log(maxfreq + EarQ*minBW) + ...
                log(minfreq + EarQ*minBW))/nfilts) * (maxfreq + EarQ*minBW);
cfreqs = flipud(cfreqs);

GTord = 4;

ucirc = exp(j*2*pi*[0:(nfft/2)]/nfft);

justpoles = 0;

for i = 1:nfilts
  cf = cfreqs(i);
  ERB = width*((cf/EarQ).^order + minBW^order).^(1/order);
  B = 1.019*2*pi*ERB;
  r = exp(-B/sr);
  theta = 2*pi*cf/sr;
  pole = r*exp(j*theta);

  if justpoles == 1
    % point on unit circle of maximum gain, from differentiating magnitude
    cosomegamax = (1+r*r)/(2*r)*cos(theta);
    if abs(cosomegamax) > 1
      if theta < pi/2;  omegamax = 0; 
      else              omegamax = pi;   end
    else
      omegamax = acos(cosomegamax);
    end
    center = exp(j*omegamax);
    gain = abs((pole-center).*(pole'-center)).^GTord;
    wts(i,1:(nfft/2+1)) = gain * (abs((pole-ucirc).*(pole'- ...
                                                     ucirc)).^-GTord);
  else
    % poles and zeros, following Malcolm's MakeERBFilter
    T = 1/sr;
    A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2* ...
                                                      cf*pi*T)./exp(B*T))/2; 
    A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2* ...
                                                      cf*pi*T)./exp(B*T))/2;
    A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2* ...
                                                      cf*pi*T)./exp(B*T))/2; 
    A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2* ...
                                                      cf*pi*T)./exp(B*T))/2; 
    zros = -[A11 A12 A13 A14]/T;
    
    gain(i) =  abs((-2*exp(4*j*cf*pi*T)*T + ...
                2*exp(-(B*T) + 2*j*cf*pi*T).*T.* ...
                (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
                 sin(2*cf*pi*T))) .* ...
               (-2*exp(4*j*cf*pi*T)*T + ...
                2*exp(-(B*T) + 2*j*cf*pi*T).*T.* ...
                (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
                 sin(2*cf*pi*T))).* ...
               (-2*exp(4*j*cf*pi*T)*T + ...
                2*exp(-(B*T) + 2*j*cf*pi*T).*T.* ...
                (cos(2*cf*pi*T) - ...
                 sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
               (-2*exp(4*j*cf*pi*T)*T + 2*exp(-(B*T) + 2*j*cf*pi*T).*T.* ...
                (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
               (-2 ./ exp(2*B*T) - 2*exp(4*j*cf*pi*T) +  ...
                2*(1 + exp(4*j*cf*pi*T))./exp(B*T)).^4);
    wts(i,1:(nfft/2+1)) = ((T^4)/gain(i)) ...
        * abs(ucirc-zros(1)).*abs(ucirc-zros(2))...
        .*abs(ucirc-zros(3)).*abs(ucirc-zros(4))...
        .*(abs((pole-ucirc).*(pole'-ucirc)).^-GTord);
  end
end

wts = wts(:,1:maxlen);