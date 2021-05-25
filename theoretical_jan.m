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

%% sats to compare
sata = [1 1];
satb = [2 1];

%% modulation
N = 1e4; %Number of bits sent
fc = 23e9;

%% constellation
P = 24; % nr of orbital planes
Q = 48; % nr of satellites per plane
h = 600e3; % altitude above the earths surface [m]
dh = 0e3; % difference in altitude between planes [m]
a = (RE+h)+dh*(1:P); % semi-major axis [m]
e = 0; % eccentricity
w = 2*pi/Q*(0:Q-1)'; % argument of periapsis w [rad]
incl = deg2rad(86.4); % inclination i [rad]
omega = pi/P*(0:P-1)'; % longitude of ascending node Omega [rad]
M0 = pi/Q*(0:P-1)'; % mean anomaly at t=0 [rad]

%% simulation steps
T = 2*pi*sqrt(a(end)^3/mu); % calculate the period of the longest orbit
Tmin = 0;
Tmax = T/2; % simulation time
dt = 60; % simulation resolution in seconds
t = (Tmin:dt:Tmax+dt)'; % make time array
L = length(t); % number of simulation steps

%% link budget variables
dEIRP = 4; %dBW/MHz- EIRP density
B = 400e6; %bandwidth [Hz]
EIRP = dEIRP+10*log10(B/1e6)+30; %dBm - EIRP
Gt = 38.5; %dBi - transmitter gain
Gr = Gt; %dBi - receiver gain
Pt = EIRP-Gt; %dBm - transmitted power
kb = 10*log10(1.3806e-23)+30; %dBm/K/Hz - boltzmann constant
Ts = 10*log10(750); %dBK - System noise temperature
Rs = 0.1*B; %Bd - expected symbol rate
Rs = 10*log10(Rs); %dBHz - symbol rate
margin = 3; %dB - other noise factors

%% estimate SNR
D = zeros(L,1);
SNR = zeros(L,1);
for k=1:L
    [ra,va] = kep2cart(a(sata(1)),e,w(sata(2)),omega(sata(1)),incl,M0(sata(1)),t(k));
    [rb,vb] = kep2cart(a(satb(1)),e,w(satb(2)),omega(satb(1)),incl,M0(satb(1)),t(k));
    D(k) = norm(ra-rb); %distance
    Lp = (4*pi*D(k)*fc/c).^2;
    Lp = 10*log10(Lp); %dB - Free space path loss
    SNR(k) = EIRP + Gr - Ts - kb - Lp - Rs - margin; %SNR per symbol
end


fig = figure();
    hold on
    box on
    grid on
    set(gca,"FontSize",14)
    xlim([min(t) max(t)]/60)
%     ylim([min(SNR)-1 max(SNR)+1])
    ylim([9 16])
    plot(t/60,SNR,'LineWidth',2)
    xlabel('Time [min]')
    ylabel('SNR per symbol [dB]')
    title("Symbol rate: "+10^(Rs/10)/1e3+" kBd")
    
%     figname="figure/linkMargin.pdf";
%     exportgraphics(fig,figname,'ContentType','vector');
%     system("pdfcrop -margins 10" + " " + figname + " " + figname);
    
    
    