function [u, x_hat] = vectorDecoder(y,const)
%% Decodes symbols into binary array using the chosen modulation
%
% Inputs: 
%   y - Arrays of symbols [complex array]
%   const - Modulation to use [char array or string]
%
% Outputs:
%   u - Decoded data [binary array]
%   x_hat - Estimated symbols [complex array]
% 
% Supported modulations:
% 	BPSK, QPSK, 8PSK, 8QAM, 16QAM, 64QAM, 256QAM
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