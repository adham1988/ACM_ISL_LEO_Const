close all
clear all
clc

%% cfg
save = true;

%% Delta
load("results/delta_closeEncounter.mat");
%time 
Tmin = 45*60; %simulation start time
dt = 10; %simulation resolution in seconds
t = Tmin:dt:Tmin+(height(results)-1)*dt;
% throughput
disp("Total throughput Delta: ")
disp(sum(results(:,8)*dt)+" Mbit")
figname1 = "results_delta.pdf";
f1 = figure();
    hold on
    box on 
    grid on
    xlim([min(t) max(t)]/60)
    xlabel('Time [min]')
    yyaxis left
    set(gca,'FontSize',14)
    ylim([0 350])
    stairs(t/60,results(:,8),'LineWidth',1.5)
    ylabel('Throughput [Mbps]')
    yyaxis right
    set(gca,'FontSize',14)
    ylim([0 1])
    stairs(t/60,results(:,4),'--','LineWidth',1.5)
    yline(1e-1,'--','LineWidth',1.5)
    ylabel('Frame Error rate')
    if save
        exportgraphics(f1,figname1,'ContentType','vector');
        system("pdfcrop -margins 10" + " " + figname1 + " " + figname1);
    end
%% Intraplane
load("results/star_intraplane.mat");
%time 
Tmin = 0; %simulation start time
dt = 10; %simulation resolution in seconds
t = Tmin:dt:Tmin+(height(results)-1)*dt;
% throughput
disp("Total throughput Intra plane: ")
disp(sum(results(:,8)*dt)+" Mbit")
figname2 = "results_intra.pdf";
f2 = figure();
    hold on
    box on 
    grid on
    set(gca,'FontSize',14)
    xlim([min(t) max(t)]/60)
    xlabel('Time [min]')
    yyaxis left
    ylim([0 350])
    stairs(t/60,results(:,8),'LineWidth',1.5)
    ylabel('Throughput [Mbps]')
    yyaxis right
    ylim([0 1])
    stairs(t/60,results(:,4),'--','LineWidth',1.5)
    yline(1e-1,'--','LineWidth',1.5)
    ylabel('Frame Error rate')
    if save
        exportgraphics(f2,figname2,'ContentType','vector');
        system("pdfcrop -margins 10" + " " + figname2 + " " + figname2);
    end
%% Interplane
load("results/star_interplane2.mat");
%time 
Tmin = 0; %simulation start time
dt = 10; %simulation resolution in seconds
t = Tmin:dt:Tmin+(height(results)-1)*dt;
% throughput
disp("Total throughput Inter plane: ")
disp(sum(results(:,8)*dt)+" Mbit")
figname3 = "results_inter.pdf";
f3 = figure();
    hold on
    box on 
    grid on
    set(gca,'FontSize',14)
    xlim([min(t) max(t)]/60)
    xlabel('Time [min]')
    yyaxis left
    ylim([0 350])
    stairs(t/60,results(:,8),'LineWidth',1.5)
    ylabel('Throughput [Mbps]')
    yyaxis right
    ylim([0 1])
    stairs(t/60,results(:,4),'--','LineWidth',1.5)
    yline(1e-1,'--','LineWidth',1.5)
    ylabel('Frame Error rate')
    if save
        exportgraphics(f3,figname3,'ContentType','vector');
        system("pdfcrop -margins 10" + " " + figname3 + " " + figname3);
    end
%% Crossplane
load("results/star_crossplane.mat");
%time 
Tmin = 0; %simulation start time
dt = 10; %simulation resolution in seconds
t = Tmin:dt:Tmin+(height(results)-1)*dt;
% throughput
disp("Total throughput cross plane: ")
disp(sum(results(:,8)*dt)+" Mbit")
figname4 = "results_cross.pdf";
f4 = figure();
    hold on
    box on 
    grid on
    set(gca,'FontSize',14)
    xlim([min(t) max(t)]/60)
    xlabel('Time [min]')
    yyaxis left
    ylim([0 350])
    stairs(t/60,results(:,8),'LineWidth',1.5)
    ylabel('Throughput [Mbps]')
    yyaxis right
    ylim([0 1])
    stairs(t/60,results(:,4),'--','LineWidth',1.5)
    yline(1e-1,'--','LineWidth',1.5)
    ylabel('Frame Error rate')
    if save
        exportgraphics(f4,figname4,'ContentType','vector');
        system("pdfcrop -margins 10" + " " + figname4 + " " + figname4);
    end









