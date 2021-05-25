more off
close all
clear all
clc
%% constants
RE = 6371e3; % radius of earth [m]
c = 299792458; % speed of light [m/s]
G = 6.67259E-11; %Gravitational constant
M = 5.9736E24; %Mass of the earth
mu = G*M; %Standard gravitational parameter
fc = 2.5e9;

%% sats to compare
sata = [1 1];
satb = [2 1];

%% message
msg_sent = "hej, det er mig!"
msg_sent = convertStringsToChars(msg_sent);
msg_sent = reshape(dec2bin(msg_sent,8).',1,[])-'0';

%% modulation
N = length(msg_sent); %Number of symbols sent
BPSK = [1; -1]; %symbol constellation
QPSK = [1 1; -1 1; -1 -1; 1 -1]'/sqrt(2); %symbol constellation
M_BPSK = length(BPSK); %number of symbols for bpsk
M_QPSK = length(QPSK); %number of symbols qpsk
K_BPSK = log2(M_BPSK); %bits per symbol
K_QPSK = log2(M_QPSK); %bits per symbol
Es = 1; %average energy per symbol

%% constellation
P = 2; % nr of orbital planes
Q = 1; % nr of satellites per plane
h = 1000e3; % altitude above the earths surface [m]
dh = 0e3; % difference in altitude between planes [m]
a = (RE+h)+dh*(1:P); % semi-major axis [m]
e = 0; % eccentricity
w = [0 0]; % argument of periapsis w [rad]
w0 = [0 0]; % argument of periapsis offset [rad]
incl = deg2rad(52); % inclination i [rad]
omega = [0 pi/8]; % longitude of ascending node Omega [rad]
M0 = 0; % mean anomaly at t=0 [rad]

%% simulation steps
P = 2*pi*sqrt(a(end)^3/mu); % calculate the period of the longest orbit
Tmin = 0;
Tmax = P/2; % simulation time
dt = 60; % simulation resolution in seconds
t = Tmin:dt:Tmax+dt; % make time array
L = length(t); % number of simulation steps

%% link budget variables
dEIRP = 34; %dBW/MHz- EIRP density
B = 30e6; %bandwidth [Hz]
EIRP = dEIRP+10*log10(B/1e6)+30; %dBm - EIRP
Gt = 30; %dBi - transmitter gain
Gr = Gt; %dBi - receiver gain
Pt = EIRP-Gt; %dBm - transmitted power
kb = 10*log10(1.3806e-23)+30; %dBm/K/Hz - boltzmann constant
Ts = 10*log10(1250); %dBK - System noise temperature
R = 10*log10(100e6); %dBHz - data rate
fudge = 2; %dB other noise factors
margin = 3+18; %dB - Signal to Noise for lossless transmission

%% THE for-loop
tic
fs = zeros(L,1);
fb = zeros(L,1);
for k=1:L
    [ra,va] = kep2cart(a(sata(1)),e,w(sata(2))+w0(sata(1)),omega(sata(1)),incl,M0,t(k));
    [rb,vb] = kep2cart(a(satb(1)),e,w(satb(2))+w0(satb(1)),omega(satb(1)),incl,M0,t(k));
    D(k) = norm(ra-rb);
    Lp = (4*pi*D(k)*fc/c).^2;
    Lp = 10*log10(Lp); %dB - Free space path loss
    SNR(k) = EIRP + Gr - Ts - kb - Lp - R - margin - fudge; % link margin
    SNR(k) = 10^(SNR(k)/10); %don't forget stupid
    N0(k) = Es/SNR(k); %noise variance
    i = 1;
    while (i<N+1)
        if(10*log10(SNR(k))<7.3)
            u = msg_sent(i);
            x = encBPSK(BPSK,u); %vector encoder
            W = sqrt(N0(k)/2)*randn(K_BPSK,1); %noise
            y = x + W; %add noise
            x_hat = invBPSK(BPSK,y); %ML decision rule
            u_hat(k,i) = decBPSK(BPSK,x_hat);
            i = i+1;
        else
            u = msg_sent([i i+1]);
            x = encQPSK(QPSK,u); %vector encoder
            W = sqrt(N0(k)/2)*randn(K_QPSK,1); %noise
            y = x + W; %add noise
            x_hat = invQPSK(QPSK,y); %ML decision rule
            logical = x==x_hat;
            u_hat(k,[i i+1]) = decQPSK(QPSK,x_hat);
            i = i+2;
        end 
    end
end

msg_recv = char(bin2dec(reshape(char(u_hat(26,:)+'0'), 8,[]).'));
msg_recv = convertCharsToStrings(msg_recv)



%% Functions
function x = encBPSK(X,u)
% BPSK encoding function
if(u==0)
    x = X(1);
elseif(u==1)
    x = X(2);
end
end

function u = decBPSK(X,x)
% BPSK decoding function
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
% QPSK decoding function
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
