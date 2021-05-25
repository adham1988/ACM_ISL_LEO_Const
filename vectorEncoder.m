function sym = vectorEncoder(u,const)
%% Encodes an array of binaries into sybols of a chosen modulation constellation
% 
% Inputs: 
%   u - Data to encode [binary array]
%   const - Modulation to use [char array or string]
%
% Outputs: 
%   sym - Array of symbols [complex array]
%
% Supported modulations:
%   BPSK, QPSK, 8PSK, 8QAM, 16QAM, 64QAM, 256QAM
sym = [];
if(const=="BPSK")
    assert(mod(length(u),1)==0,'Expected u to be divisible by 1 for BPSK')
    for i=1:length(u)
        data = bi2de(u(i));
        sym = [sym pskmod(data,2,0)];
    end
elseif(const=="QPSK")
    assert(mod(length(u),2)==0,'Expected u to be divisible by 2 for QPSK')
    for i=1:2:length(u)
        data = bi2de(reshape(u(i:i+1),1,[]));
        sym = [sym pskmod(data,4,pi/4)];
    end
elseif(const=="8PSK")
    assert(mod(length(u),3)==0,'Expected u to be divisible by 3 for 8PSK')
    for i=1:3:length(u)
        data = bi2de(reshape(u(i:i+2),1,[]));
        sym = [sym pskmod(data,8,pi/8)];
    end
elseif(const=="8QAM")
    assert(mod(length(u),3)==0,'Expected u to be divisible by 3 for 8QAM')
    for i=1:3:length(u)
        data = bi2de(reshape(u(i:i+2),1,[]));
        sym = [sym qammod(data,8,'unitAveragePower',true)];
    end
elseif(const=="16QAM")
    assert(mod(length(u),4)==0,'Expected u to be divisible by 4 for 16QAM')
    for i=1:4:length(u)
        data = bi2de(reshape(u(i:i+3),1,[]));
        sym = [sym qammod(data,16,'unitAveragePower',true)];
    end
elseif(const=="64QAM")
    assert(mod(length(u),6)==0,'Expected u to be divisible by 6 for 64QAM')
    for i=1:6:length(u)
        data = bi2de(reshape(u(i:i+5),1,[]));
        sym = [sym qammod(data,64,'unitAveragePower',true)];
    end
elseif(const=="256QAM")
    assert(mod(length(u),8)==0,'Expected u to be divisible by 8 for 256QAM')
    for i=1:8:length(u)
        data = bi2de(reshape(u(i:i+7),1,[]));
        sym = [sym qammod(data,256,'unitAveragePower',true)];
    end
end
sym = reshape(sym,[],1);
end