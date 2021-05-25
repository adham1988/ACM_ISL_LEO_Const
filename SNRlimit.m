more off
close all
clear all
clc

%% Constants
fc = 23e9;
B = 400e6;
Rs = 0.1*B;
%% Colors
colors = [  0,      0.4470, 0.741;
            0.8500, 0.3250, 0.0980;
            0.9290, 0.6940, 0.1250;
            0.4940, 0.1840, 0.5560;
            0.4660, 0.6740, 0.1880;
            0.3010, 0.7450, 0.9330;
            0.6350, 0.0780, 0.1840;]';
%% INIT
SNR = 0:0.01:40; %SNR per symbol
K = [1 2 3 4 6 8]; %bits per symbol
p = 3:9; %number of parity bits
k = 2.^(p)-p-1; %number of data bits
n = p+k; %codeword length
Rc = k./n; %code rate
L = 1536; %frame length

%% Uncoded
Uncoded = zeros(length(SNR),length(K));
Uncoded(:,1) = berawgn(SNR-10*log10(K(1)),'psk',2^K(1),'nondiff');
Uncoded(:,2) = berawgn(SNR-10*log10(K(2)),'psk',2^K(2),'nondiff');
Uncoded(:,3) = berawgn(SNR-10*log10(K(3)),'psk',2^K(3),'nondiff');
for i=4:length(K)
    Uncoded(:,i) = berawgn(SNR-10*log10(K(i)),'qam',2^K(i));
end

%% Hamming
Hamming = zeros(length(SNR),length(K),length(p));
for j=1:length(n)
    Hamming(:,1,j) = bercoding(SNR-10*log10(K(1)),'hamming','hard',n(j),'psk',2^K(1),'nondiff');
    Hamming(:,2,j) = bercoding(SNR-10*log10(K(2)),'hamming','hard',n(j),'psk',2^K(2),'nondiff');
    Hamming(:,3,j) = bercoding(SNR-10*log10(K(3)),'hamming','hard',n(j),'psk',2^K(3),'nondiff');
    for i=4:length(K)
        Hamming(:,i,j) = bercoding(SNR-10*log10(K(i)),'hamming','hard',n(j),'qam',2^K(i));
    end
end

%% Bit error Probability
BER = zeros(length(SNR),length(K),length(p));
BER(:,:,1:7) = Hamming;
BER(:,:,end+1) = Uncoded;

fig1 = figure();
    hold on
    box on
    grid on
    set(gca,'yscale','log')
    set(gca,'fontsize',14)
    ylim([1e-8 1e0])
    xlim([min(SNR) max(SNR)])
    for i=1:length(K)
        for j=1:length(p)+1
            plot(SNR,BER(:,i,j),'color',colors(:,i))
        end
    end
    yline(1e-5,'k--')
    xlabel('SNR per symbol [dB]')
    ylabel('BER')
    
%% Frame error probability
FER = 1-(1-(BER)).^L;
figname2 = "figure/PER.pdf";
fig2 = figure();
    hold on
    box on
    grid on
    set(gca,'yscale','log')
    set(gca,'fontsize',14)
    ylim([1e-8 1e0])
    xlim([min(SNR) max(SNR)])
    for i=1:length(K)
        for j=1:length(p)+1
            plot(SNR,FER(:,i,j),'color',colors(:,i))
        end
    end
    yline(1e-2,'k--')
    xlabel('SNR per symbol [dB]')
    ylabel('FER')
    exportgraphics(fig2,figname2,'ContentType','vector');
    system("pdfcrop -margins 10" + " " + figname2 + " " + figname2);
    
%% limits
lim = 1e-2;
SNR_req = zeros(length(K),length(p)+1);
for i=1:length(K)
    for j=1:length(p)+1
        SNR_lim = SNR(FER(:,i,j)<lim);
        SNR_lim = SNR_lim(1);
        SNR_req(i,j) = SNR_lim;
    end
end
p = [p 0];
Rc = [Rc 1];
%%
figname3 = "figure/stair.pdf";
fig3 = figure();
    hold on
    box on
    grid on
    xlim([8 9])
    ylim([15 45])
    set(gca,'FontSize',14)
    for j=1:length(p)
        stairs([0 SNR_req(:,j)'],[0 Rs.*K(:)'.*Rc(j)/1e6],'-','LineWidth',1.5)
    end
    for i=1:length(K)
        for j=1:length(p)
            text(SNR_req(i,j),Rs*K(i)*Rc(j)/1e6+1,"p:"+p(j))
        end
    end
    xlabel('Required SNR per symbol [dB]')
    ylabel('Throughput [Mbps]')
    title("FER: "+lim)
    exportgraphics(fig3,figname3,'ContentType','vector');
    system("pdfcrop -margins 10" + " " + figname3 + " " + figname3);

%% SNR lookup table
Throughput = zeros(length(K),length(p));
for i=1:length(K)
    for j=1:length(p)
        Throughput(i,j) = Rs*K(i)*Rc(j)/1e6;
    end
end
I = [1];
J = [4];
Throughput_skift = Throughput(I,J);
SNR_skift = SNR_req(I,J);
for m=1:length(SNR)
    for i=1:length(K)
        for j=1:length(p)
            if(SNR(m)>SNR_req(i,j))
                if(Throughput(i,j)>Throughput(I(end),J(end)))
                    I = [I i];
                    J = [J j];
                    Throughput_skift = [Throughput_skift Throughput(i,j)];
                    SNR_skift = [SNR_skift SNR_req(i,j)];
                end
            end
        end
    end
end

lookupTable = [K(I)' p(J)' SNR_skift' Throughput_skift']

for i=1:length(I)
    if(K(I(i))==1)
        label(i) = "BPSK";
    elseif(K(I(i))==2)
        label(i) = "QPSK";
    elseif(K(I(i))==3)
        label(i) = "8QAM";
    elseif(K(I(i))==4)
        label(i) = "16QAM";
    elseif(K(I(i))==6)
        label(i) = "64QAM";
    elseif(K(I(i))==8)
        label(i) = "256QAM";
    end
end
%%
filename="lookUptable.csv";
fid = fopen(filename,'w');
fprintf(fid,"Mod,K,p,SNR,Throughput\n");
for i=1:length(I)
    fprintf(fid,'%s,%i,%i,%f,%f\n',label(i),K(I(i)),p(J(i)),SNR_skift(i),Throughput_skift(i));
end
fclose(fid);
