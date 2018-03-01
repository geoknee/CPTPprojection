% read in ensemble_size running times and precisions for each method, for each d.
clear;close all;
dmax = 5;
ensemble_size = 20;
gdapB_times  = zeros(ensemble_size,dmax);
mosek_times  = zeros(ensemble_size,dmax);
sdpt3_times  = zeros(ensemble_size,dmax);
DIA_times    = zeros(ensemble_size,dmax);
% sedumi_times = zeros(ensemble_size,dmax); % this fails

gdapB_errors    = zeros(ensemble_size,dmax);
mosek_errors    = zeros(ensemble_size,dmax);
sdpt3_errors    = zeros(ensemble_size,dmax);
DIA_errors      = zeros(ensemble_size,dmax);
% sedumi_errors   = zeros(ensemble_size,dmax);


for d=2:dmax
    for i=1:ensemble_size
        
        dir = sprintf('./benchmarking_results/d%i',d);
        clear choi_ground
        load([dir,'/dataset',num2str(i)]);
        
        clear choi_ml_vec
        
        load([dir,'/DIA_results',num2str(i)]);
        DIA_times(i,d)  = elapsedTime;
        choi_DIA        = reshape(choi_ml_vec,[],d*d);
        DIA_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_DIA/trace(choi_DIA));
        
        clear choi_ml_vec
               
        load([dir,'/gdapB_results',num2str(i)]);
        gdapB_times(i,d) = elapsedTime;
        choi_gdapB        = reshape(choi_ml_vec,[],d*d);
        gdapB_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_gdapB/trace(choi_gdapB));
        
        clear choi_ml_vec
        
        load([dir,'/mosek_results',num2str(i)]);
        mosek_times(i,d) = elapsedTime;    
        choi_mosek        = reshape(choi_ml_vec,[],d*d);
        mosek_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_mosek/trace(choi_mosek));
        
        clear choi_ml_vec
        
        load([dir,'/sdpt3_results',num2str(i)]);
        sdpt3_times(i,d) = elapsedTime;    
        choi_sdpt3        = reshape(choi_ml_vec,[],d*d);
        sdpt3_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_sdpt3/trace(choi_sdpt3));
        
%         clear choi_ml_vec
%         
%         load([dir,'/sedumi_results',num2str(i)]);
%         sedumi_times(i,d) = elapsedTime;    
%         choi_sedumi        = reshape(choi_ml_vec,[],d*d);
%         sedumi_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_sedumi/trace(choi_sedumi));
        
        
    end
end

figure; hold on;

errorbar(2:dmax,mean(DIA_times(:,2:end)),std(DIA_times(:,2:end)),'-d','LineWidth',2);
errorbar(2:dmax,mean(gdapB_times(:,2:end)),std(gdapB_times(:,2:end)),'-*','LineWidth',2);
errorbar(2:dmax,mean(mosek_times(:,2:end)),std(mosek_times(:,2:end)),'-x','LineWidth',2);
errorbar(2:dmax,mean(sdpt3_times(:,2:end)),std(sdpt3_times(:,2:end)),'-x','LineWidth',2);
% errorbar(2:dmax,mean(sedumi_times(:,2:end)),std(sedumi_times(:,2:end)),'-s','LineWidth',2);
xlim([1.8,dmax+0.2])
xlabel 'Hilbert space dimension'
ylabel 'times taken (s)';
set(gca,'YScale','log');
legend('DIA','gdapB','mosek','sdpt3','sedumi')
% legend('gdapB','mosek','sdpt3')
box on
grid on
set(gca,'fontsize',18)
saveas(gcf,'time.png')
saveas(gcf,'time.eps','epsc')


figure; hold on;

errorbar(2:dmax,mean(DIA_errors(:,2:end)),std(DIA_errors(:,2:end)),'-d','LineWidth',2);
errorbar(2:dmax,mean(gdapB_errors(:,2:end)),std(gdapB_errors(:,2:end)),'-*','LineWidth',2);
errorbar(2:dmax,mean(mosek_errors(:,2:end)),std(mosek_errors(:,2:end)),'-x','LineWidth',2);
errorbar(2:dmax,mean(sdpt3_errors(:,2:end)),std(sdpt3_errors(:,2:end)),'-x','LineWidth',2);
% errorbar(2:dmax,mean(sedumi_errors(:,2:end)),std(sedumi_errors(:,2:end)),'-s','LineWidth',2);
xlim([1.8,dmax+0.2])
ylim([0,1])
xlabel 'Hilbert space dimension'
ylabel 'J distance';
set(gca,'YScale','log')
legend('DIA','gdapB','mosek','sdpt3','sedumi')
% legend('gdapB','mosek','sdpt3')
box on
grid on
set(gca,'fontsize',18)
saveas(gcf,'errors.png')
saveas(gcf,'errors.eps','epsc')



figure; hold on;
scatter(DIA_errors(:,2),DIA_times(:,2),'o','filled','DisplayName','DIA')
scatter(gdapB_errors(:,2),gdapB_times(:,2),'o','filled','DisplayName','gdapB')
scatter(mosek_errors(:,2),mosek_times(:,2),'o','filled','DisplayName','mosek')
scatter(sdpt3_errors(:,2),sdpt3_times(:,2),'o','filled','DisplayName','sdpt3')
% scatter(sedumi_errors(:,2),sedumi_times(:,2),'o','filled','DisplayName','sedumi')

legend('show','Location','northwest')

ax = gca;
ax.ColorOrderIndex = 1;

scatter(DIA_errors(:,3),DIA_times(:,3),'<','filled')
scatter(gdapB_errors(:,3),gdapB_times(:,3),'<','filled')
scatter(mosek_errors(:,3),mosek_times(:,3),'<','filled')
scatter(sdpt3_errors(:,3),sdpt3_times(:,3),'<','filled')
% scatter(sedumi_errors(:,3),sedumi_times(:,3),'<','filled')

ax = gca;
ax.ColorOrderIndex = 1;

scatter(DIA_errors(:,4),DIA_times(:,4),'d','filled')
scatter(gdapB_errors(:,4),gdapB_times(:,4),'d','filled')
scatter(mosek_errors(:,4),mosek_times(:,4),'d','filled')
scatter(sdpt3_errors(:,4),sdpt3_times(:,4),'d','filled')
% scatter(sedumi_errors(:,4),sedumi_times(:,4),'d','filled')

ax = gca;
ax.ColorOrderIndex = 1;

scatter(DIA_errors(:,5),DIA_times(:,5),'p','filled')
scatter(gdapB_errors(:,5),gdapB_times(:,5),'p','filled')
scatter(mosek_errors(:,5),mosek_times(:,5),'p','filled')
scatter(sdpt3_errors(:,5),sdpt3_times(:,5),'p','filled')
% scatter(sedumi_errors(:,5),sedumi_times(:,5),'p','filled')

ax = gca;
ax.ColorOrderIndex = 1;

scatter(DIA_errors(:,6),DIA_times(:,6),'s','filled')
scatter(gdapB_errors(:,6),gdapB_times(:,6),'s','filled')
scatter(mosek_errors(:,6),mosek_times(:,6),'s','filled')
scatter(sdpt3_errors(:,6),sdpt3_times(:,6),'s','filled')
% scatter(sedumi_errors(:,6),sedumi_times(:,6),'s','filled')

set(gca,'xscale','log')
set(gca,'yscale','log')

xlabel 'error'
ylabel 'time (s)'

box on
grid on



set(gca,'fontsize',18)
saveas(gcf,'scatter.png')
saveas(gcf,'scatter.eps','epsc')
