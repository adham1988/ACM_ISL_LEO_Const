%load("results/delta_closeEncounter.mat");
%load("results/star_intraplane.mat");
%load("results/star_interplane2.mat");
load("results/star_crossplane.mat");

newResults = results;

for b = 1:length(newResults)
    Hlen = 11;
    p = results(b,7);
    k = (2^p)-1-p;
    n = (2^p);
    N = 1536/n;
    R = ((k*N)-Hlen)/(N*n);
    if p == 0
        R = (1536-11)/1536;
    end
    
    Rs = 40;
    
    newResults(b,10) = Rs*newResults(b,6)*R*(1-newResults(b,4));
    
end

%time
Tmin = 0; %simulation start time
dt = 10; %simulation resolution in seconds
t = Tmin:dt:Tmin+(height(newResults)-1)*dt;
% throughput
disp("Total throughput cross plane: ")
disp(sum(newResults(:,8)*dt)+" Mbit")
figname1 = "results_cross_goodput.pdf";
f1 = figure();
hold on
box on
grid on
xlim([min(t) max(t)]/60)
xlabel('Time [min]')
yyaxis right
set(gca,'FontSize',14)
ylim([0 350])
stairs(t/60,newResults(:,8),'--','LineWidth',1.5)
ylabel('Throughput [Mbps]')
yyaxis left
set(gca,'FontSize',14)
ylim([0 350])
stairs(t/60,newResults(:,10),'LineWidth',1.5)
ylabel('Goodput [Mbps]')
    exportgraphics(f1,figname1,'ContentType','vector');
    system("pdfcrop -margins 10" + " " + figname1 + " " + figname1);
