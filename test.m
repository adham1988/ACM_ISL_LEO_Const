close all
clc

%%

M = 64;
k = log2(M);
q = quantizer([k 0],'ufixed');

data = randi([0 1],1,k)
data = bi2de(data);
sym = qammod(data,M);
data = qamdemod(sym,M);
data = dec2bin(data,k);
nyData = [];
for i=1:length(data)
    nyData(length(data)-i+1) = str2double(data(i));
end

nyData