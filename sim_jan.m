more off
close all
clear all
clc

%% constants
RE = 6371e3; %radius of earth [m]
c = 299792458; %speed of light [m/s]
G = 6.67259E-11; %Gravitational constant
M = 5.9736E24; %Mass of the earth
mu = G*M; %Standard gravitational parameter
fc = 23e9; %carrier frequency

%% sats to compare
sata = [1 1];
satb = [2 1];

%% message
% txtFile = "message.txt";
% FID = fopen(txtFile);
% msg_sent = textscan(FID,'%s','delimiter','');
% fclose(FID);
% msg_sent = string(msg_sent{:});
% for i=2:length(msg_sent)
%     msg_sent(1) = msg_sent(1) + msg_sent(i) + "\n";
% end
% msg_sent = msg_sent(1);
% % disp(sprintf("message sent: \n\n"+msg_sent+"\n"))
% msg_sent = dec2bin(char(msg_sent),8);
% msg_sent = msg_sent(:,end-7:end);
% msg_sent = reshape(msg_sent.',1,[])-'0';
msg_sent = double(rand(1536*10,1)<0.5)';

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
dt = 60; %simulation resolution in seconds
t = Tmin:dt:Tmax+dt; %make time array
L = length(t); %number of simulation steps

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

%% modulation
N = length(msg_sent); %Number of bits sent
%BPSK
BPSK = [1; -1]; %symbol constellation
%QPSK
QPSK = [cos(1*pi/4) sin(1*pi/4); cos(3*pi/4) sin(3*pi/4); ...
        cos(5*pi/4) sin(5*pi/4); cos(7*pi/4) sin(7*pi/4);]'; 
%8PSK
EightPSK = [cos(1*pi/4) sin(1*pi/4); cos(2*pi/4) sin(2*pi/4); ...
            cos(3*pi/4) sin(3*pi/4); cos(4*pi/4) sin(4*pi/4); ...
            cos(5*pi/4) sin(5*pi/4); cos(6*pi/4) sin(6*pi/4); ...
            cos(7*pi/4) sin(7*pi/4); cos(8*pi/4) sin(8*pi/4);]';
%8QAM
EightQAM = [1+sqrt(3) 0; 1 1; 1 -1; 0 1+sqrt(3); ...
            -1 1; -1 -1; 0 -1-sqrt(3); -1-sqrt(3) 0;]';
s = 0;
for i=1:length(EightQAM)
    s = s + norm(EightQAM(:,i));
end
s = s/length(EightQAM);
EightQAM = EightQAM/s;
%16QAM
SixteenQAM = [-3 3; -1 3; 1 3; 3 3; ...
              -3 1; -1 1; 1 1; 3 1; ...
              -3 -1; -1 -1; 1 -1; 3 -1; ...
              -3 -3; -1 -3; 1 -3; 3 -3;]';
SixteenQAM = SixteenQAM/sqrt(10);

%% make higher order QAM constellations
SixtyfourQAM = makeConst(64)/sqrt(42);
TwohundredfiftysixQAM = makeConst(256)/sqrt(170);

%% modulation constants
M_BPSK = length(BPSK); %number of symbols for bpsk
M_QPSK = length(QPSK); %number of symbols qpsk
M_EightPSK = length(EightPSK); %number of symbols 8psk
M_EightQAM = length(EightQAM); %number of symbols 8QAM
M_SixteenQAM = length(SixteenQAM); %number of symbols 16QAM
M_SixtyfourQAM = length(SixtyfourQAM); %number of symbols 64QAM
M_TwohundredfiftysixQAM = length(TwohundredfiftysixQAM); %number of symbols 256QAM

K_BPSK = log2(M_BPSK); %bits per symbol
K_QPSK = log2(M_QPSK); %bits per symbol
K_EightPSK = log2(M_EightPSK); %bits per symbol
K_EightQAM = log2(M_EightQAM); %bits per symbol
K_SixteenQAM = log2(M_SixteenQAM); %bits per symbol
K_SixtyfourQAM = log2(M_SixtyfourQAM); %bits per symbol
K_TwohundredfiftysixQAM = log2(M_TwohundredfiftysixQAM); %bits per symbol

%% THE for-loop
tic
D = zeros(L,1);
SNR = zeros(L,1);
N0 = zeros(L,1);
u_hat = zeros(L,N);
Pb = zeros(L,1);
parfor k=1:L
    tic
    [ra,va] = kep2cart(a(sata(1)),e,w(sata(2)),omega(sata(1)),incl,M0(sata(1)),t(k));
    [rb,vb] = kep2cart(a(satb(1)),e,w(satb(2)),omega(satb(1)),incl,M0(satb(1)),t(k));
    D(k) = norm(ra-rb);
    Lp = (4*pi*D(k)*fc/c).^2;
    Lp = 10*log10(Lp); %dB - Free space path loss
    SNR(k) = EIRP + Gr - Ts - kb - Lp - R - margin; %link margin
    SNR(k) = 10^(SNR(k)/10); %don't forget stupid
    N0(k) = Es/SNR(k); %noise variance
    i = 1;
    N0_est = [];
    
    
    if(10*log10(SNR(k))>30.48)
        const = "256QAM";
    elseif(10*log10(SNR(k))>24.54)
        const = "64QAM";
    elseif(10*log10(SNR(k))>18.64)
        const = "16QAM";
    elseif(10*log10(SNR(k))>16.30)
        const = "8PSK";
    elseif(10*log10(SNR(k))>11.63)
        const = "QPSK";
    elseif(10*log10(SNR(k))>8.60)
        const = "BPSK";
    end
    
    u = msg_sent;
    x = vectorEncoder(u,const);
    y = channel(x,N0(k)*ones(length(x),1));
    [u_hat(k,:),x_hat] = vectorDecoder(y,const);
    N0_est = channelEstimator(y,x_hat);
    SNR_est(k) = Es/mean(N0_est);
    Pb(k) = sum(u_hat(k,:)~=msg_sent,2)/N;
    msg_recv(k,:) = char(bin2dec(reshape(char(u_hat(k,:)+'0'), 8,[]).'));
    toc
end
toc

%% Error probabilities
%BPSK
Ps_BPSK = 1/2*erfc(sqrt(SNR));
Pb_BPSK = Ps_BPSK;
%QPSK
Ps_QPSK = erfc(sqrt(SNR));
Pb_QPSK = 1/2*erfc(sqrt(SNR/2));
%8PSK
Ps_EightPSK = erfc(sqrt(SNR)*sin(pi/M_EightPSK));
Pb_EightPSK = Ps_EightPSK/K_EightPSK;
%8QAM;
Ps_EightQAM = 2*(1-1/sqrt(M_EightQAM))*erfc(sqrt(3/(2*(M_EightQAM-1))*SNR));
Pb_EightQAM = Ps_EightQAM/K_EightQAM;
%16QAM;
Ps_SixteenQAM = 2*(1-1/sqrt(M_SixteenQAM))*erfc(sqrt(3/(2*(M_SixteenQAM-1))*SNR));
Pb_SixteenQAM = Ps_SixteenQAM/K_SixteenQAM;
%64QAM;
Ps_SixtyfourQAM = 2*(1-1/sqrt(M_SixtyfourQAM))*erfc(sqrt(3/(2*(M_SixtyfourQAM-1))*SNR));
Pb_SixtyfourQAM = Ps_SixtyfourQAM/K_SixtyfourQAM;
%16QAM;
Ps_TwohundredfiftysixQAM = 2*(1-1/sqrt(M_TwohundredfiftysixQAM))*erfc(sqrt(3/(2*(M_TwohundredfiftysixQAM-1))*SNR));
Pb_TwohundredfiftysixQAM = Ps_TwohundredfiftysixQAM/K_TwohundredfiftysixQAM;

%% received message
msg_recv = convertCharsToStrings(msg_recv(1,:));
% disp(sprintf("message received: \n\n"+msg_recv+"\n"))

%% plot
% figure
%     hold on
%     box on
%     grid on
%     set(gca,'yscale','log')
%     xlim([min(t) max(t)]/60)
%     plot(t/60,Ps,'*')
%     plot(t/60,Ps_BPSK,'-')
%     plot(t/60,Ps_QPSK,'-')
%     plot(t/60,Ps_EightPSK,'-')
%     xlabel('Time [min]')
%     ylabel('Ps')
%     legend('Measured','BPSK','QPSK','8PSK')

figure
    hold on
    box on
    grid on
    set(gca,'yscale','log')
    xlim([min(t) max(t)]/60)
    ylim([1e-6 1e0])
    plot(t/60,Pb,'*')
    plot(t/60,Pb_BPSK,'-')
    plot(t/60,Pb_QPSK,'-')
    plot(t/60,Pb_EightPSK,'-')
    plot(t/60,Pb_SixteenQAM,'-')
    plot(t/60,Pb_SixtyfourQAM,'-')
    plot(t/60,Pb_TwohundredfiftysixQAM,'-')
    xlabel('Time [min]')
    ylabel('Pb')
    legend('Measured','BPSK','QPSK','8PSK','16QAM','64QAM','256QAM')
    
figure
    hold on
    box on
    grid on
    set(gca,'yscale','log')
    ylim([1e-6 1e0])
    plot(10*log10(SNR),Pb,'*')
    plot(10*log10(SNR),Pb_BPSK,'-')
    plot(10*log10(SNR),Pb_QPSK,'-')
    plot(10*log10(SNR),Pb_EightPSK,'-')
    plot(10*log10(SNR),Pb_SixteenQAM,'-')
    plot(10*log10(SNR),Pb_SixtyfourQAM,'-')
    plot(10*log10(SNR),Pb_TwohundredfiftysixQAM,'-')
    xlabel('SNR [dB]')
    ylabel('Pb')
    legend('Measured','BPSK','QPSK','8PSK','16QAM','64QAM','256QAM')
    
figure
    hold on
    box on
    grid on
    plot(t/60,10*log10(SNR),'-')
    plot(t/60,10*log10(SNR_est),'*')
    xlabel('Time [min]')
    ylabel('SNR per symbol [dB]')
    legend('SNR','Estimated SNR')
    
    
%% Functions
function x = encBPSK(X,u)
%BPSK encoding function
if(u==0)
    x = X(1);
elseif(u==1)
    x = X(2);
end
end
function u = decBPSK(X,x)
%BPSK decoding function
if(x==X(1))
    u = 0;
elseif(x==X(2))
    u = 1;
end
end
function x = invBPSK(X,y)
%BPSK encoding function
if(y>=0)
    x = X(1);
elseif(y<0)
    x = X(2);
end
end
function x = encQPSK(X,u)
%QPSK encoding function
if (u(1)==0 && u(2)==0)
    x = X(:,1);
elseif (u(1)==0 && u(2)==1)
    x = X(:,2);
elseif (u(1)==1 && u(2)==1)
    x = X(:,3);
elseif (u(1)==1 && u(2)==0)
    x = X(:,4);
end
end
function u = decQPSK(X,x)
%QPSK decoding function
if(x(1)==X(1,1) && x(2)==X(2,1))
    u = [0 0]';
elseif(x(1)==X(1,2) && x(2)==X(2,2))
    u = [0 1]';
elseif(x(1)==X(1,3) && x(2)==X(2,3))
    u = [1 1]';
elseif(x(1)==X(1,4) && x(2)==X(2,4))
    u = [1 0]';
end
end
function x = invQPSK(X,y)
%QPSK ML decision rule
if(y(1)>=0 && y(2)>=0)
    x = X(:,1);
elseif(y(1)<0 && y(2)>=0)
    x = X(:,2);
elseif(y(1)<0 && y(2)<0)
    x = X(:,3);
elseif(y(1)>=0 && y(2)<0)
    x = X(:,4);
end

% M = zeros(length(X),1);
% I = zeros(length(X),1);
% for i=1:length(X)
%     [M(i),I(i)] = min(norm(X(:,i)-y)); %decision rule
% end
% [m,j] = min(M);
% x = X(:,j); %choose symbol
end
function x = encEightPSK(X,u)
%QPSK encoding function
if(u(1)==0 && u(2)==0 && u(3)==0)
    x = X(:,1);
elseif(u(1)==0 && u(2)==0 && u(3)==1)
    x = X(:,2);
elseif(u(1)==0 && u(2)==1 && u(3)==1)
    x = X(:,3);
elseif(u(1)==0 && u(2)==1 && u(3)==0)
    x = X(:,4);
elseif(u(1)==1 && u(2)==1 && u(3)==0)
    x = X(:,5);
elseif(u(1)==1 && u(2)==1 && u(3)==1)
    x = X(:,6);
elseif(u(1)==1 && u(2)==0 && u(3)==1)
    x = X(:,7);
elseif(u(1)==1 && u(2)==0 && u(3)==0)
    x = X(:,8);
end
end
function u = decEightPSK(X,x)
%QPSK decoding function
if(x(1)==X(1,1) && x(2)==X(2,1))
    u = [0 0 0]';
elseif(x(1)==X(1,2) && x(2)==X(2,2))
    u = [0 0 1]';
elseif(x(1)==X(1,3) && x(2)==X(2,3))
    u = [0 1 1]';
elseif(x(1)==X(1,4) && x(2)==X(2,4))
    u = [0 1 0]';
elseif(x(1)==X(1,5) && x(2)==X(2,5))
    u = [1 1 0]';
elseif(x(1)==X(1,6) && x(2)==X(2,6))
    u = [1 1 1]';
elseif(x(1)==X(1,7) && x(2)==X(2,7))
    u = [1 0 1]';
elseif(x(1)==X(1,8) && x(2)==X(2,8))
    u = [1 0 0]';
end
end
function x = invQAM(X,y)
%QAM ML decision rule
M = zeros(length(X),1);
I = zeros(length(X),1);
for i=1:length(X)
    [M(i),I(i)] = min(norm(X(:,i)-y)); %decision rule
end
[m,j] = min(M);
x = X(:,j); %choose symbol
end
function x = encEightQAM(X,u)
%8QAM encoding function
if(u(1)==0 && u(2)==0 && u(3)==0)
    x = X(:,1);
elseif(u(1)==0 && u(2)==0 && u(3)==1)
    x = X(:,2);
elseif(u(1)==1 && u(2)==0 && u(3)==0)
    x = X(:,3);
elseif(u(1)==0 && u(2)==1 && u(3)==1)
    x = X(:,4);
elseif(u(1)==0 && u(2)==1 && u(3)==0)
    x = X(:,5);
elseif(u(1)==1 && u(2)==1 && u(3)==1)
    x = X(:,6);
elseif(u(1)==1 && u(2)==0 && u(3)==1)
    x = X(:,7);
elseif(u(1)==1 && u(2)==1 && u(3)==0)
    x = X(:,8);
end
end
function u = decEightQAM(X,x)
%8QAM decoding function
if(x(1)==X(1,1) && x(2)==X(2,1))
    u = [0 0 0]';
elseif(x(1)==X(1,2) && x(2)==X(2,2))
    u = [0 0 1]';
elseif(x(1)==X(1,3) && x(2)==X(2,3))
    u = [1 0 0]';
elseif(x(1)==X(1,4) && x(2)==X(2,4))
    u = [0 1 1]';
elseif(x(1)==X(1,5) && x(2)==X(2,5))
    u = [0 1 0]';
elseif(x(1)==X(1,6) && x(2)==X(2,6))
    u = [1 1 1]';
elseif(x(1)==X(1,7) && x(2)==X(2,7))
    u = [1 0 1]';
elseif(x(1)==X(1,8) && x(2)==X(2,8))
    u = [1 1 0]';
end
end
function x = encSixteenQAM(X,u)
%16QAM encoding function
if(u(1)==0 && u(2)==0 && u(3)==0 && u(4)==0)
    x = X(:,1);
elseif(u(1)==0 && u(2)==1 && u(3)==0 && u(4)==0)
    x = X(:,2);
elseif(u(1)==1 && u(2)==1 && u(3)==0 && u(4)==0)
    x = X(:,3);
elseif(u(1)==1 && u(2)==0 && u(3)==0 && u(4)==0)
    x = X(:,4);
elseif(u(1)==0 && u(2)==0 && u(3)==0 && u(4)==1)
    x = X(:,5);
elseif(u(1)==0 && u(2)==1 && u(3)==0 && u(4)==1)
    x = X(:,6);
elseif(u(1)==1 && u(2)==1 && u(3)==0 && u(4)==1)
    x = X(:,7);
elseif(u(1)==1 && u(2)==0 && u(3)==0 && u(4)==1)
    x = X(:,8);
elseif(u(1)==0 && u(2)==0 && u(3)==1 && u(4)==1)
    x = X(:,9);
elseif(u(1)==0 && u(2)==1 && u(3)==1 && u(4)==1)
    x = X(:,10);
elseif(u(1)==1 && u(2)==1 && u(3)==1 && u(4)==1)
    x = X(:,11);
elseif(u(1)==1 && u(2)==0 && u(3)==1 && u(4)==1)
    x = X(:,12);
elseif(u(1)==0 && u(2)==0 && u(3)==1 && u(4)==0)
    x = X(:,13);
elseif(u(1)==0 && u(2)==1 && u(3)==1 && u(4)==0)
    x = X(:,14);
elseif(u(1)==1 && u(2)==1 && u(3)==1 && u(4)==0)
    x = X(:,15);
elseif(u(1)==1 && u(2)==0 && u(3)==1 && u(4)==0)
    x = X(:,16);
end
end
function u = decSixteenQAM(X,x)
%16QAM decoding function
if(x(1)==X(1,1) && x(2)==X(2,1))
    u = [0 0 0 0]';
elseif(x(1)==X(1,2) && x(2)==X(2,2))
    u = [0 1 0 0]';
elseif(x(1)==X(1,3) && x(2)==X(2,3))
    u = [1 1 0 0]';
elseif(x(1)==X(1,4) && x(2)==X(2,4))
    u = [1 0 0 0]';
elseif(x(1)==X(1,5) && x(2)==X(2,5))
    u = [0 0 0 1]';
elseif(x(1)==X(1,6) && x(2)==X(2,6))
    u = [0 1 0 1]';
elseif(x(1)==X(1,7) && x(2)==X(2,7))
    u = [1 1 0 1]';
elseif(x(1)==X(1,8) && x(2)==X(2,8))
    u = [1 0 0 1]';
elseif(x(1)==X(1,9) && x(2)==X(2,9))
    u = [0 0 1 1]';
elseif(x(1)==X(1,10) && x(2)==X(2,10))
    u = [0 1 1 1]';
elseif(x(1)==X(1,11) && x(2)==X(2,11))
    u = [1 1 1 1]';
elseif(x(1)==X(1,12) && x(2)==X(2,12))
    u = [1 0 1 1]';
elseif(x(1)==X(1,13) && x(2)==X(2,13))
    u = [0 0 1 0]';
elseif(x(1)==X(1,14) && x(2)==X(2,14))
    u = [0 1 1 0]';
elseif(x(1)==X(1,15) && x(2)==X(2,15))
    u = [1 1 1 0]';
elseif(x(1)==X(1,16) && x(2)==X(2,16))
    u = [1 0 1 0]';
end
end
function Const = makeConst(M)
L = sqrt(M);
for j = 1:L
    for i = 1:L
        a(j,i) = -L + 1 + 2*(i-1);
        b(i,j) = -L + 1 + 2*(i-1);
    end
end
for i = 1:M
   c(:,i) = [b(i) a(i)]; 
end
Const = c;
end

%% physics
function [D,V] = dist(r,v,a,b)
% calculates the distance and radial velocity between two satellites (a and b)
%
% Inputs:
% r,v - position and velocity matrix with indicies (time,plane,nr,coord)
% a,b - satellite index vectors [plane, nr]
%
% Outputs: 
% D - distance between a and b at all times specified in r
% V - radial velocity between a and b at all times specified in r
D = zeros(length(r(:,1,1,1)),1);
V = zeros(length(v(:,1,1,1)),1);
for i=1:length(r(:,1,1,1))
    D(i) = norm(squeeze(r(i,b(1),b(2),:))-squeeze(r(i,a(1),a(2),:)));
    V(i) = norm((squeeze(v(i,b(1),b(2),:))-squeeze(v(i,a(1),a(2),:)))' ...
        *(squeeze(r(i,b(1),b(2),:))-squeeze(r(i,a(1),a(2),:))))/D(i);
end
end
function [r,v] = kep2cart(a,e,w,omega,i,M0,t)
% Calculates Cartesian state vectors from a set of Kelperian orbit elements
% 
% Inputs: 
% a - semi major axis [m]
% e - eccentricity
% w - argument of periapsis [rad]
% omega - longitude of ascending node [rad]
% i - inclination [rad]
% M0 - mean anomaly at time 0
% t - time
% 
% Outputs:
% r - cartesian position vector at time t [m]
% v - cartesian velocity vector at time t [m/s]
G = 6.67259E-11; %Gravitational constant
M = 5.9736E24; %Mass of the earth
mu = G*M; %Standard gravitational parameter

n = sqrt(mu/a^3); %rate of sweep
M = mod(M0 + n*t,2*pi); %mean anomaly at time t

E = solveKeplerEq(M,e,1e-10); %solve kepler equation for eccentric anomaly

nu = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2)); %calculate true anomaly from eccentric anomaly

rc = a*(1-e*cos(E)); %distance from center of gravity

o = rc * [cos(nu) sin(nu) 0]'; %position in kepler coordinates
odot = sqrt(mu*a)/rc * [-sin(E) sqrt(1-e^2)*cos(E) 0]'; %velocity in kepler coordinates

%transform position to cartesian reference frame
r = [o(1)*(cos(w)*cos(omega)-sin(w)*cos(i)*sin(omega)) - o(2)*(sin(w)*cos(omega)+cos(w)*cos(i)*sin(omega));
    o(1)*(cos(w)*sin(omega)+sin(w)*cos(i)*cos(omega)) + o(2)*(cos(w)*cos(i)*cos(omega)-sin(w)*sin(omega));
    o(1)*(sin(w)*sin(i)) + o(2)*(cos(w)*sin(i))];
%transform velocity to cartesian reference frame
v = [odot(1)*(cos(w)*cos(omega)-sin(w)*cos(i)*sin(omega)) - odot(2)*(sin(w)*cos(omega)+cos(w)*cos(i)*sin(omega));
    odot(1)*(cos(w)*sin(omega)+sin(w)*cos(i)*cos(omega)) + odot(2)*(cos(w)*cos(i)*cos(omega)-sin(w)*sin(omega));
    odot(1)*(sin(w)*sin(i)) + odot(2)*(cos(w)*sin(i))];
end
function E = solveKeplerEq(M,e,tol)
% Solves Kepler's equation M(t) = E(t) - e*sin(E(t)) for E(t) via the
% Newton-Raphson method. 
% 
% Inputs: 
% M - Mean anomaly at time t [rad]
% e - eccentricity
% tol - desired numerical accuracy
%
% Outputs: 
% E - eccentric anomaly at time t [rad]
E = M;
E(end+1)= E(end) - (E(end)-e*sin(E(end))-M)/(1-e*cos(E(end)));
while( abs(E(end)-E(end-1)) > tol )
    E(end+1) = E(end) - (E(end)-e*sin(E(end))-M)/(1-e*cos(E(end)));
end
E = E(end);
end