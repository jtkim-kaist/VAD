function y = sigmoid(y, fac)
% SIGMOID nonlinear funcion for cochlear model
%	y = sigmoid(y, fac);
%	fac: non-linear factor
%	 -- fac > 0, transister-like function
%	 -- fac = 0, hard-limiter
%	 -- fac = -1, half-wave rectifier
%	 -- else, no operation, i.e., linear 
%
%	SIGMOID is a monotonic increasing function which simulates 
%	hair cell nonlinearity. 
%	See also: WAV2AUD, AUD2WAV
 
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

if fac > 0,
	%y = exp(y/fac); y = 1/(1+.1)-1./(1+.1*y);
	y = exp(-y/fac); y = 1./(1+y);
elseif fac == 0,
	y = (y > 0);	% hard-limiter
elseif fac == -1,
	y = max(y, 0);	% half-wave rectifier
elseif fac == -3,
	y = halfregu(y);
end;
