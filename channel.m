function y = channel(x,N0)
%% Adds white gaussian noise (awgn) to an array of symbols
%
% Inputs: 
%   x - Symbols [complex array]
%   N0 - Noise spectral density [fraction] (NOT dB)
% 
% Outputs: 
%   y - Symbols with added noise
x = reshape(x,[],1);
N0 = reshape(N0,[],1);
assert(length(x)==length(N0),'x must be the same length as N0')
sigma = sqrt(N0/2).*randn(length(N0),2);
W = [];
for i=1:length(N0)
    W(i) = complex(sigma(i,1),sigma(i,2));
end
W = reshape(W,[],1);
y = x + W;
end