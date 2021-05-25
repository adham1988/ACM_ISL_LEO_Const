more off
close all
clear all
clc

%%
data = randi([0 1],1,1536);
const = "BPSK";

sym = vectorEncoder(data,const);

SNR = 9
N0 = 0*(1:length(sym))'+1/10^(SNR/10);
y = channel(sym,N0);

scatterplot(y)
grid

[u,x_hat] = vectorDecoder(y,const);

SNR_est = 10*log10(1/channelEstimator(y,x_hat))

Errors = sum(data~=u);
Ps = Errors/length(data);



%% Channel estimator
function N0 = channelEstimator(y,x)
y = reshape(y,[],1);
x = reshape(x,[],1);
N0 = [];
assert(length(y)==length(x),'y and x must be the same length')
for i=1:length(y)
    N0 = [2*abs(real(y)-real(x)).^2 2*abs(imag(y)-imag(x)).^2];
end
N0 = mean(reshape(N0,[],1));
end
%% Decoder
function [u, x_hat] = vectorDecoder(y,const)
u = [];
x_hat = [];
y = reshape(y,[],1);
if(const=="BPSK")
    for i=1:length(y)
        data = pskdemod(y(i),2,0);
        data = dec2bin(data,1);
        nyData = [];
        for j=1:length(data)
            nyData(length(data)-j+1) = str2double(data(j));
        end
        x_hat(i) = pskmod(bi2de(nyData),2,0);
        u = [u nyData];
    end
elseif(const=="QPSK")
    for i=1:length(y)
        data = pskdemod(y(i),4,pi/4);
        data = dec2bin(data,2);
        nyData = [];
        for j=1:length(data)
            nyData(length(data)-j+1) = str2double(data(j));
        end
        x_hat(i) = pskmod(bi2de(nyData),4,pi/4);
        u = [u nyData];
    end
elseif(const=="8PSK")
    for i=1:length(y)
        data = pskdemod(y(i),8,pi/8);
        data = dec2bin(data,3);
        nyData = [];
        for j=1:length(data)
            nyData(length(data)-j+1) = str2double(data(j));
        end
        x_hat(i) = pskmod(bi2de(nyData),8,pi/8);
        u = [u nyData];
    end
elseif(const=="8QAM")
    for i=1:length(y)
        data = qamdemod(y(i),8,'unitAveragePower',true);
        data = dec2bin(data,3);
        nyData = [];
        for j=1:length(data)
            nyData(length(data)-j+1) = str2double(data(j));
        end
        x_hat(i) = qammod(bi2de(nyData),8,'unitAveragePower',true);
        u = [u nyData];
    end
elseif(const=="16QAM")
    for i=1:length(y)
        data = qamdemod(y(i),16,'unitAveragePower',true);
        data = dec2bin(data,4);
        nyData = [];
        for j=1:length(data)
            nyData(length(data)-j+1) = str2double(data(j));
        end
        x_hat(i) = qammod(bi2de(nyData),16,'unitAveragePower',true);
        u = [u nyData];
    end
elseif(const=="64QAM")
    for i=1:length(y)
        data = qamdemod(y(i),64,'unitAveragePower',true);
        data = dec2bin(data,6);
        nyData = [];
        for j=1:length(data)
            nyData(length(data)-j+1) = str2double(data(j));
        end
        x_hat(i) = qammod(bi2de(nyData),64,'unitAveragePower',true);
        u = [u nyData];
    end
elseif(const=="256QAM")
    for i=1:length(y)
        data = qamdemod(y(i),256,'unitAveragePower',true);
        data = dec2bin(data,8);
        nyData = [];
        for j=1:length(data)
            nyData(length(data)-j+1) = str2double(data(j));
        end
        x_hat(i) = qammod(bi2de(nyData),256,'unitAveragePower',true);
        u = [u nyData];
    end
end
x_hat = reshape(x_hat,[],1);
end
%% Channel
function y = channel(x,N0)
x = reshape(x,[],1);
N0 = reshape(N0,[],1);
assert(length(x)==length(N0),'x must be the same length as N0')
sigma = sqrt(N0/2).*randn(length(N0),2);
W = [];
for i=1:length(sigma)
    W(i,:) = complex(sigma(i,1),sigma(i,2));
end
y = x + W;
end
%% Encoder
function x = vectorEncoder(u,const)
x = [];
if(const=="BPSK")
    assert(mod(length(u),1)==0,'Expected u to be divisible by 1')
    for i=1:length(u)
        data = bi2de(u(i));
        x = [x pskmod(data,2,0)];
    end
elseif(const=="QPSK")
    assert(mod(length(u),2)==0,'Expected u to be divisible by 2')
    for i=1:2:length(u)
        data = bi2de(reshape(u(i:i+1),1,[]));
        x = [x pskmod(data,4,pi/4)];
    end
elseif(const=="8PSK")
    assert(mod(length(u),3)==0,'Expected u to be divisible by 3')
    for i=1:3:length(u)
        data = bi2de(reshape(u(i:i+2),1,[]));
        x = [x pskmod(data,8,pi/8)];
    end
elseif(const=="8QAM")
    assert(mod(length(u),3)==0,'Expected u to be divisible by 3')
    for i=1:3:length(u)
        data = bi2de(reshape(u(i:i+2),1,[]));
        x = [x qammod(data,8,'unitAveragePower',true)];
    end
elseif(const=="16QAM")
    assert(mod(length(u),4)==0,'Expected u to be divisible by 4')
    for i=1:4:length(u)
        data = bi2de(reshape(u(i:i+3),1,[]));
        x = [x qammod(data,16,'unitAveragePower',true)];
    end
elseif(const=="64QAM")
    assert(mod(length(u),6)==0,'Expected u to be divisible by 6')
    for i=1:6:length(u)
        data = bi2de(reshape(u(i:i+5),1,[]));
        x = [x qammod(data,64,'unitAveragePower',true)];
    end
elseif(const=="256QAM")
    assert(mod(length(u),8)==0,'Expected u to be divisible by 8')
    for i=1:8:length(u)
        data = bi2de(reshape(u(i:i+7),1,[]));
        x = [x qammod(data,256,'unitAveragePower',true)];
    end
end
x = reshape(x,[],1);
end