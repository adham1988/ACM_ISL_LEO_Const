more off
close all
clear all
clc

%% 
SNR = 0:0.01:30;
SNR = 10.^(SNR/10)';
L = 1536;
K = [1 2 3 4 6 8 10];
CR = [1 4/7 247/255 233/255];
fc = 23e9;
B = 400e6;
Rs = 0.1*B;

%% uncoded
U_BPSK = berawgn(10*log10(SNR),'psk',2,'nondiff');
U_QPSK = berawgn(10*log10(SNR/2),'psk',4,'nondiff');
U_8PSK = berawgn(10*log10(SNR/3),'psk',8,'nondiff');
U_16QAM = berawgn(10*log10(SNR/4),'qam',16);
U_64QAM = berawgn(10*log10(SNR/6),'qam',64);
U_256QAM = berawgn(10*log10(SNR/8),'qam',256);
U_1024QAM = berawgn(10*log10(SNR/10),'qam',1024);

Uncoded = [U_BPSK U_QPSK U_8PSK U_16QAM U_64QAM U_256QAM U_1024QAM];
PEP_Uncoded = zeros(length(Uncoded),width(Uncoded));
for i=1:width(Uncoded)
    PEP_Uncoded(:,i) = 1-(1-Uncoded(:,i)).^L;
end
G_Uncoded = zeros(length(Uncoded),width(Uncoded));
for i=1:width(Uncoded)
    G_Uncoded(:,i) = K(i)*Rs.*(1-PEP_Uncoded(:,i));
end
%% Hamming 7
n = 7;
H7_BPSK = bercoding(10*log10(SNR),'hamming','hard',n,'psk',2,'nondiff');
H7_QPSK = bercoding(10*log10(SNR/2),'hamming','hard',n,'psk',4,'nondiff');
H7_8PSK = bercoding(10*log10(SNR/3),'hamming','hard',n,'psk',8,'nondiff');
H7_16QAM = bercoding(10*log10(SNR/4),'hamming','hard',n,'qam',16);
H7_64QAM = bercoding(10*log10(SNR/6),'hamming','hard',n,'qam',64);
H7_256QAM = bercoding(10*log10(SNR/8),'hamming','hard',n,'qam',256);
H7_1024QAM = bercoding(10*log10(SNR/10),'hamming','hard',n,'qam',1024);

Hamming7 = [H7_BPSK H7_QPSK H7_8PSK H7_16QAM H7_64QAM H7_256QAM H7_1024QAM];
PEP_Hamming7 = zeros(length(Hamming7),width(Hamming7));
for i=1:width(Hamming7)
    PEP_Hamming7(:,i) = 1-(1-Hamming7(:,i)).^L;
end
G_Hamming7 = zeros(length(Hamming7),width(Hamming7));
for i=1:width(Hamming7)
    G_Hamming7(:,i) = K(i)*Rs.*(1-PEP_Hamming7(:,i))*4/7;
end
%% Hamming 255
n = 255;
H255_BPSK = bercoding(10*log10(SNR),'hamming','hard',n,'psk',2,'nondiff');
H255_QPSK = bercoding(10*log10(SNR/2),'hamming','hard',n,'psk',4,'nondiff');
H255_8PSK = bercoding(10*log10(SNR/3),'hamming','hard',n,'psk',8,'nondiff');
H255_16QAM = bercoding(10*log10(SNR/4),'hamming','hard',n,'qam',16);
H255_64QAM = bercoding(10*log10(SNR/6),'hamming','hard',n,'qam',64);
H255_256QAM = bercoding(10*log10(SNR/8),'hamming','hard',n,'qam',256);
H255_1024QAM = bercoding(10*log10(SNR/10),'hamming','hard',n,'qam',1024);

Hamming255 = [H255_BPSK H255_QPSK H255_8PSK H255_16QAM H255_64QAM H255_256QAM H255_1024QAM];
PEP_Hamming255 = zeros(length(Hamming255),width(Hamming255));
for i=1:width(Hamming255)
    PEP_Hamming255(:,i) = 1-(1-Hamming255(:,i)).^L;
end
G_Hamming255 = zeros(length(Hamming255),width(Hamming255));
for i=1:width(Hamming255)
    G_Hamming255(:,i) = K(i)*Rs.*(1-PEP_Hamming255(:,i))*247/255;
end
%% Reed-solomon
% n = 255;
% k = 233;
% RS_BPSK = bercoding(10*log10(SNR),'RS','hard',n,k,'psk',2,'nondiff');
% RS_QPSK = bercoding(10*log10(SNR/2),'RS','hard',n,k,'psk',4,'nondiff');
% RS_8PSK = bercoding(10*log10(SNR/3),'RS','hard',n,k,'psk',8,'nondiff');
% RS_16QAM = bercoding(10*log10(SNR/4),'RS','hard',n,k,'qam',16);
% RS_64QAM = bercoding(10*log10(SNR/6),'RS','hard',n,k,'qam',64);
% RS_256QAM = bercoding(10*log10(SNR/8),'RS','hard',n,k,'qam',256);
% RS_1024QAM = bercoding(10*log10(SNR/10),'RS','hard',n,k,'qam',1024);
% 
% ReedSolomon = [RS_BPSK RS_QPSK RS_8PSK RS_16QAM RS_64QAM RS_256QAM RS_1024QAM];
% PEP_ReedSolomon = zeros(length(ReedSolomon),width(ReedSolomon));
% for i=1:width(ReedSolomon)
%     PEP_ReedSolomon(:,i) = 1-(1-ReedSolomon(:,i)).^L;
% end
% G_ReedSolomon = zeros(length(ReedSolomon),width(ReedSolomon));
% for i=1:width(ReedSolomon)
%     G_ReedSolomon(:,i) = K(i)*Rs.*(1-PEP_ReedSolomon(:,i))*233/255;
% end

%% Colors
colors = [0, 0.4470, 0.741;
            0.8500, 0.3250, 0.0980;
            0.9290, 0.6940, 0.1250;
            0.4940, 0.1840, 0.5560;
            0.4660, 0.6740, 0.1880;
            0.3010, 0.7450, 0.9330;
            0.6350, 0.0780, 0.1840;]';
        
%% BER
figure
    hold on
    box on
    grid on
    set(gca,'Yscale','log')
    ylim([1e-10 1e0])
    %uncoded
    for i=1:width(Uncoded)
        plot(10*log10(SNR),Uncoded(:,i),'-','color',colors(:,i),'LineWidth',2)
    end
    %hamming7
    for i=1:width(Hamming7)
        plot(10*log10(SNR),Hamming7(:,i),'--','color',colors(:,i),'LineWidth',2)
    end
    %hamming255
    for i=1:width(Hamming255)
        plot(10*log10(SNR),Hamming255(:,i),':','color',colors(:,i),'LineWidth',2)
    end
%     %reed-solomon
%     for i=1:width(ReedSolomon)
%         plot(10*log10(SNR),ReedSolomon(:,i),'-.','color',colors(:,i),'LineWidth',2)
%     end
    yline(1e-5,'k--')
    xlabel('SNR [dB]')
    ylabel('BER')
    legend('BPSK','QPSK','8PSK','16QAM','64QAM','256QAM','1024QAM')

%% PEP
figure
    hold on
    box on
    grid on
    set(gca,'Yscale','log')
    ylim([1e-3 1e0])
    %uncoded
    for i=1:width(Uncoded)
        plot(10*log10(SNR),PEP_Uncoded(:,i),'-','color',colors(:,i),'LineWidth',2)
    end
    %hamming7
    for i=1:width(Hamming7)
        plot(10*log10(SNR),PEP_Hamming7(:,i),'--','color',colors(:,i),'LineWidth',2)
    end
    %hamming255
    for i=1:width(Hamming255)
        plot(10*log10(SNR),PEP_Hamming255(:,i),':','color',colors(:,i),'LineWidth',2)
    end
%     %reed-solomon
%     for i=1:width(ReedSolomon)
%         plot(10*log10(SNR),PEP_ReedSolomon(:,i),'-.','color',colors(:,i),'LineWidth',2)
%     end
    yline(1e-1,'k--')
    xlabel('SNR [dB]')
    ylabel('PER')
    legend('BPSK','QPSK','8PSK','16QAM','64QAM','256QAM','1024QAM')
    
    
    
    
    
    return
%% Throughput

lim = 1e-1;
SNR_req = zeros(width(Uncoded),length(CR));
for i=1:width(Uncoded)
    U_lim = PEP_Uncoded(:,i);
    SNR_lim = SNR(U_lim<lim);
    SNR_lim = SNR_lim(1);
    SNR_req(i,1) = SNR_lim;
end
for i=1:width(Hamming7)
    U_lim = PEP_Hamming7(:,i);
    SNR_lim = SNR(U_lim<lim);
    SNR_lim = SNR_lim(1);
    SNR_req(i,2) = SNR_lim;
end
for i=1:width(Hamming255)
    U_lim = PEP_Hamming255(:,i);
    SNR_lim = SNR(U_lim<lim);
    SNR_lim = SNR_lim(1);
    SNR_req(i,3) = SNR_lim;
end
for i=1:width(ReedSolomon)
    U_lim = PEP_ReedSolomon(:,i);
    SNR_lim = SNR(U_lim<lim);
    SNR_lim = SNR_lim(1);
    SNR_req(i,4) = SNR_lim;
end

figure
    hold on
    box on
    grid on
    ylim([0 200])
    for i=1:3
        stairs([0 10*log10(SNR_req(:,i))'],[0 (Rs*K*CR(i))/1e6])
    end
    xlabel('Required SNR per symbol [dB]')
    ylabel('Throughput [Mbps]')
    legend('Uncoded','H7','H255','RS')
    title("PER: "+lim)














