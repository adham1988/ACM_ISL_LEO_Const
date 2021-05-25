function N0 = channelEstimator(y,x)
%% Estimates the noise spectral density from the received and estimated symbols
% 
% Inputs:
%   y - Received symbols [complex array]
%   x - estimated symbols [complex array]
%
% Outputs:
%   N0 - Estimated noise spectral density [fraction] (NOT dB)
y = reshape(y,[],1);
x = reshape(x,[],1);
N0 = [];
assert(length(y)==length(x),'y and x must be the same length')
for i=1:length(y)
    N0 = [2*abs(real(y)-real(x)).^2 2*abs(imag(y)-imag(x)).^2];
end
N0 = mean(reshape(N0,[],1));
end