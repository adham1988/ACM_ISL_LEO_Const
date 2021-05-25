clear all
clc


Pt  = 1;            %set power transmittet in dBW
Gr  = 30;           %set receiver gain in dB
Gt  = 30;           %set transmitter gain in dB
fc  = 23e9;         %center frequency Hz
B   = 100e6;        %bandwidth Hz
c   = 299792458;    %speed of light m/s
Kb  = 10*log10(1.3806e-23);     %boltzmann constant dBW/K/Hz
Ts  = 10*log10(1250);           %Noise temperature dBK
Rs  = 10*log10([0.01 0.005 0.0015 0.001 0.00025 0.000062 0.000016]*B);            %symbolrate in dBHz
j   = [1 2 3 4 6 8 10]; %bits pr symbol
d   = 1000e3;                   %fixed distance

sata = [1 1];
satb = [2 1];
%% constants
RE = 6371e3; % radius of earth [m]
c = 299792458; % speed of light [m/s]
G = 6.67259E-11; %Gravitational constant
M = 5.9736E24; %Mass of the earth
mu = G*M; %Standard gravitational parameter
 % Example
    N = 2; % nr of orbital planes
	M = 1; % nr of satellites per plane
    h = 1000e3; % altitude above the earths surface [m]
    dh = 0e3; % difference in altitude between planes [m]
    a = (RE+h)+dh*(1:N); % semi-major axis [m]
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
dt = 1; % simulation resolution in seconds
t = (Tmin:dt:Tmax+dt)'; % make time array
L = length(t); % number of simulation steps
%% estimate SNR
D = zeros(L,1);
SNR = zeros(L,1);
req = [9.6  12.6 17.8 19.5  25.6  31.6   37.5]; 
     %[BPSK QPSK 8psk 16QAM 64QAM 256QAM 1024QAM]
for k=1:L
    [ra,va] = kep2cart(a(sata(1)),e,w(sata(2))+w0(sata(1)),omega(sata(1)),incl,M0,t(k));
    [rb,vb] = kep2cart(a(satb(1)),e,w(satb(2))+w0(satb(1)),omega(satb(1)),incl,M0,t(k));
    D(k) = norm(ra-rb); %distance
    Lp = (4*pi*D(k)*fc/c).^2;
    Lp = 10*log10(Lp); %dB - Free space path loss
for i = 1:length(req)
    SNR(k, i) = Pt + Gt + Gr - Ts - Kb - Lp - Rs(i) - req(i); %SNR per symbol
end
end
%%
Lp = (4*pi*d*fc/c)^2;   %Path loss
Lp = 10*log10(Lp);      %now in dB

figure
hold on
grid on
box on
plot(t, SNR)
% title("Symbolrate "+10^(Rs/10)/1e6+" MBd" )
xlabel("tid [s]")
ylabel("Margin [dB]")
xlim([min(t) max(t)])
legend("BPSK", "QPSK", "8PSK", "16QAM", "64QAM" ,"256QAM", "1024QAM")



















