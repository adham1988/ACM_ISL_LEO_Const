close all
clc

%%

%% constants
RE = 6371e3; %radius of earth [m]
c = 299792458; %speed of light [m/s]
G = 6.67259E-11; %Gravitational constant
M = 5.9736E24; %Mass of the earth
mu = G*M; %Standard gravitational parameter
fc = 23e9; %carrier frequency
tid = 0;

%% sats to compare
sata = [1 1];
satb = [2 1];

%% constellation
P = 24; %nr of orbital planes
Q = 48; %nr of satellites per plane
h = 600e3; %altitude above the earths surface [m]
dh = 0e3; %difference in altitude between planes [m]
a = (RE+h)+dh*(1:P); %semi-major axis [m]
e = 0; %eccentricity
w = 2*pi/Q*(0:Q-1)'; %argument of periapsis w [rad]
incl = deg2rad(86.4); %inclination i [rad]
omega = pi/P*(0:P-1)'; %longitude of ascending node Omega [rad]
M0 = pi/Q*(0:P-1)'; %mean anomaly at t=0 [rad]

%% simulation steps
T = 2*pi*sqrt(a(end)^3/mu); %calculate the period of the longest orbit
Tmin = 0; %simulation start time
Tmax = T/2; %simulation stop time
dt = 10; %simulation resolution in seconds
t = Tmin:dt:Tmax+dt; %make time array
L = length(t); %number of simulation steps

%% BSS generation and frame parameters
dataPerFrame(502) = 3*502;
dataPerFrame(247) = 6*247;
dataPerFrame(120) = 12*120;
dataPerFrame(57) = 24*57;
dataPerFrame(26) = 48*26;
dataPerFrame(11) = 96*11;
dataPerFrame(4) = 192*4;

bssSize = 5e4; % Amount of bytes from the Binary Source Stream.
lenHeader = zeros(1,11); % Length of length header
N = 1536; % Maximum frame size containing the header and coded bites.
frameAmount = 100; % Amount of frames sent per step
allowedChars = ['0':'1'];
whichChars = randi(numel(allowedChars), [1 bssSize]); % Chars of size P in Bytes
bss = allowedChars(whichChars); % Binary source stream
bss = convertStringsToChars(bss);
bss = reshape(dec2bin(bss,8).',1,[])-'0';
LUT = readmatrix("lookUptable.csv"); % import LUT

%% link budget variables
dEIRP = 4; %dBW/MHz- EIRP density
B = 400e6; %bandwidth [Hz]
EIRP = dEIRP+10*log10(B/1e6)+30; %dBm - EIRP
Gt = 38.5; %dBi - transmitter gain
Gr = Gt; %dBi - receiver gain
Pt = EIRP-Gt; %dBm - transmitted power
kb = 10*log10(1.3806e-23)+30; %dBm/K/Hz - boltzmann constant
Ts = 10*log10(750); %dBK - System noise temperature
R = 10*log10(0.1*B); %dBHz - data rate
margin = 3; %dB - Signal to Noise for lossless transmission
Es = 1;

%% THE for-loop
tic
D = zeros(L,1);
SNR = zeros(L,1);
SNR_est = zeros(L,1);
SNR_est(1) = 7.17;
N0 = zeros(L,1);
Pb = zeros(L,1);
BIS = double.empty;
results = zeros(L,9); % Vector of results. 0 errors uncoded | 1 errors uncoded | 2 errors uncoded | PER | SNR

for k=1:L
    [ra,va] = kep2cart(a(sata(1)),e,w(sata(2)),omega(sata(1)),incl,M0(sata(1)),t(k));
    [rb,vb] = kep2cart(a(satb(1)),e,w(satb(2)),omega(satb(1)),incl,M0(satb(1)),t(k));
    D(k) = norm(ra-rb);
    Lp = (4*pi*D(k)*fc/c).^2;
    Lp = 10*log10(Lp); %dB - Free space path loss
    SNR(k) = EIRP + Gr - Ts - kb - Lp - R - margin; %link margin
    SNR(k) = 10^(SNR(k)/10); %don't forget stupid
    N0(k) = Es/SNR(k); %noise variance
    
    if k == 1
        [symbolSize,parityLength,LUTSNR,throughput] = ACM_LUT(LUT,7.17);
    else
        % Choose parameters from Lookup Table based on SNR
        [symbolSize,parityLength,LUTSNR,throughput] = ACM_LUT(LUT,10*log10(SNR_est(k-1)));
    end

    datawordLength = 2^parityLength - 1 - parityLength; % k. Data length = 2^m-1-m
    codewordLength = datawordLength + parityLength; % n. Codeword length = 2^m-1

    u_hat = zeros(L,N);

    for g = 1:frameAmount % Repeat for the amount of packets
        i = 1; % Reset i

        % Encode a frame
        if parityLength > 0 % Encode if SNR specifies any parity
            [G,H,frame_sent,backupFrame,bss] = encoding(bss,datawordLength, ...
                lenHeader,dataPerFrame,codewordLength,parityLength);
        else % Uncoded
            if length(bss) >= 1525
                frameData = bss(1:1525);
                headerLen = dec2bin(length(frameData),11);
                frame = [headerLen,frameData];
                bss = bss(length(frameData):end);
                frame_sent = mod(double(frame),2);
                backupFrame = mod(double(frame),2);
            else
                frameData = bss(1:end);
                headerLen = dec2bin(length(frameData),11);
                frame = [headerLen frameData zeros(1,N-length(frameData)-length(headerLen))];
                bss = [];
                frame_sent = mod(double(frame),2);
                backupFrame = mod(double(frame),2);
            end
        end

        % Modulation and demodulation

        const = chooseMod(symbolSize); %ACM controller
        x = vectorEncoder(frame_sent,const); %modulation
        y = channel(x,N0(k)*ones(length(x),1)); %apply channel
        [u_hat(k,:),x_hat] = vectorDecoder(y,const); %demodulate
        N0_est = channelEstimator(y,x_hat); %estimate channel
        SNR_est(k,g) = Es/mean(N0_est);
        % Modulation and demodulation end

        Pb(k) = sum(u_hat(k,:)~=frame_sent,2)/N;
        msg_recv = char(u_hat(k,:)+'0');

        % Decoding encoded
        if datawordLength > 0 % Decodes coded
            [results,BIS] = decoding(msg_recv,H,G,datawordLength,codewordLength,backupFrame,results,k,BIS);
        else % Uncoded
            if mod(double(msg_recv),2) == mod(backupFrame,2)
                results(k,1) = results(k,1) + 1;
                lengthofFrameBits = bin2dec(msg_recv(1:11));
                BIS = [BIS msg_recv(11:lengthofFrameBits+11)];
            else
                results(k,3) = results(k,3) + 1;
            end
        end
    end

    results(k,4) = 1-(results(k,1)/frameAmount);
    if(k==1)
        results(k,5) = 7.17;
    else
        results(k,5) = 10*log10(mean(SNR_est(k,:)));
    end
    results(k,6) = symbolSize;
    results(k,7) = parityLength;
    results(k,8) = throughput;
    results(k,9) = 10*log10(SNR(k));
    disp(k)
end

header = ["None" "One" "Two" "PER" "estimatedSNR" "Mod" "Par" "Throughput" "trueSNR"];

toc
save star_intraplane.mat results

%% plots
figname1 = "estimatesSNR_star_intraplane.pdf";
f1 = figure();
    hold on
    box on
    grid on
    set(gca,'FontSize',14)
    xlim([min(t) max(t)]/60)
    plot(t/60,10*log10(SNR),'LineWidth',1)
    stairs(t/60,results(:,5),'LineWidth',1)
    xlabel('Time [min]')
    ylabel('SNR per symbol [dB]')
    legend('True','Estimated')
    exportgraphics(f1,figname1,'ContentType','vector');
    system("pdfcrop -margins 10" + " " + figname1 + " " + figname1);

%% functions 
% chose modulation from number of bits per symbol from LUT
function const = chooseMod(symsize)
if(symsize==8)
    const = "256QAM";
elseif(symsize==6)
    const = "64QAM";
elseif(symsize==4)
    const = "16QAM";
elseif(symsize==3)
    const = "8PSK";
elseif(symsize==2)
    const = "QPSK";
elseif(symsize==1)
    const = "BPSK";
else
    const = "None";
end
end

function [symbolSize,paritySize,LUTSNR,throughput] = ACM_LUT(LUT,SNR)
l = 0;
while l < height(LUT)
    l = l + 1;
    if(l==length(LUT))
        symbolSize = LUT(1,2);
        paritySize = LUT(1,3);
        LUTSNR = LUT(1,4);
        throughput = LUT(1,5);
        l = l + 50;
    else
        if SNR-LUT(end+1-l,4) > 0
            symbolSize = LUT(end+1-l,2);
            paritySize = LUT(end+1-l,3);
            LUTSNR = LUT(end+1-l,4);
            throughput = LUT(end+1-l,5);
            l = l + 50;
        end
    end
end
end

% Encoding
function [G,H,frame_sent,backupFrame,bss] = encoding(bss,datawordLength,lenHeader,dataPerFrame,codewordLength,parityLength)

frame_sent = double.empty; % Empty array for storing transmission databits
if length(bss) >= dataPerFrame(datawordLength)
    frameData = bss(1:dataPerFrame(datawordLength)-length(lenHeader));
    headerLen = dec2bin(length(frameData),11);
    frame = [headerLen frameData];
    bss = bss(length(frameData):end);
else
    frameData = bss(1:end);
    headerLen = dec2bin(length(frameData),11);
    frame = [headerLen frameData zeros(1,dataPerFrame(datawordLength)-length(frameData)-length(headerLen))];
    bss = [0];
end
backupFrame = mod(double(frame),2);
[H,G] = generateMatrices(codewordLength,datawordLength,parityLength); % Generate matrices
while ~isempty(frame) % Encode entire frame
    dataword = frame(1:datawordLength); % Fill a temporary dataword
    codeword = mod(dataword*G,2); % c. Codeword generation
    if mod(sum(codeword),2) == 1 % Extended parity
        codeword(end+1) = 1; % Even parity as a one
    else
        codeword(end+1) = 0; % Even parity as a zero
    end
    frame_sent = [frame_sent,codeword]; % Add to transmitter frame
    frame = frame(datawordLength+1:end); % Remove bits from frame
end
end

% Decoding
function [results,BIS] = decoding(msg_recv,H,G,datawordLength,codewordLength,backupFrame,results,k,BIS)
remErrors = 0; % Detected errors
msg_decoded = double.empty; % Empty array for storing receiver databits
while ~isempty(msg_recv) % Until packet is empty...
    r = mod(double(msg_recv(1:codewordLength+1)),2);
    syndrome = mod(r * H.',2); % s. Syndrome
    if sum(syndrome) == 0 % Syndrome is a null vector, meaning no errors.
        msg_decoded = [msg_decoded,r(:,1:datawordLength)]; % Add to received packet
        msg_recv = msg_recv(codewordLength+2:end); % Remove sent bits from packet
    else           % Syndrome is not a null vector, error was detected.
        [cv,cvErr] = codewordGen(datawordLength,codewordLength,G,r); % Estimate the valid codeword
        if cvErr == 0 % If correct codeword was found
            msg_decoded = [msg_decoded,cv(:,1:datawordLength)]; % Add to received packet
            msg_recv = msg_recv(codewordLength+2:end); % Remove sent bits from packet
        elseif cvErr == 1 % Two or more errors
            msg_decoded = [msg_decoded,r(:,1:datawordLength)]; % Add to received packet
            msg_recv = msg_recv(codewordLength+2:end); % Remove sent bits from packet
            remErrors = 1; % An error remains
        end
    end
end
if remErrors == 1 % Remainng error known, discard frame
    results(k,2) = results(k,2) + 1;
else
    if mod(msg_decoded,2) == mod(backupFrame,2)
        lengthofFrameBits = bin2dec(char(msg_decoded(1:11)+'0'));
        BIS = [BIS msg_decoded(11:lengthofFrameBits+11)];
        results(k,1) = results(k,1) + 1;
    else
        results(k,3) = results(k,3) + 1;
    end
end
end

% Generate Generator and Parity Check Matrices
% H = Parity Check Matrix
% G = Generator Matrix
function [H,G] = generateMatrices(n,k,m)
[~,g] = hammgen(m); % Create scaled parity matrix
p = g(:,1:m); % k x n-k parity matrix
clearvars g % Free some mem
IG = eye(k); % k x k Identity matrix for Generator Matrix
IH = eye(n-k); % n-k x n-k Identity matrix for Parity Check Matrix
G = [IG p]; % Generator Matrix
H = [p' IH]; % Parity Check Matrix
H(:, end+1) = 0;
H = [H ; ones(1,n+1)];
%z = sum(mod(G*H',2),'all'); % G*H^T should give 0.
clearvars IG IH % Free some mem
end

% Checks if single bit changes make a valid codeword.
% cv = Generated valid codeword
% cvErr = 0 for corrected error, 1 for more than one error
function [cv,cvErr] = codewordGen(k,n,G,r)
for j = 1:n
    %disp(j)
    r(j) = mod(r(j)+1,2); % Flip one bit
    cv = mod(r(:,1:k)*G,2); % cv. Generate valid code from data bits
    if mod(sum(cv),2) == 1
        cv(end+1) = 1; % Even parity as a one
    else
        cv(end+1) = 0; % Even parity as a zero
    end
    distance = biterr(r,cv); % Compare received codeword with the valid codeword
    if distance == 0 % Hamming Distance = 1 for valid codeword
        %disp("Correct codeword found");
        cvErr = 0;
        return
    end
    r(j) = mod(r(j)-1,2); % Flip it back
end
%disp("More than one error detected");
cvErr = 1;
end
