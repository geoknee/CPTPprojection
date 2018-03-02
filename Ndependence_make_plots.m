% read in ensemble_size running times and precisions for each method, for each d.
clear;close all;
dmax = 4;
ensemble_size = 10;
Npows = [1,2,3,4,5,6,7,8,9,Inf];
Ns = 10.^Npows;
Ns(end)=10^12;
% Ns = [2,4,8,16,32];
gdapB_times  = zeros(ensemble_size,dmax,length(Ns));
mosek_times  = zeros(ensemble_size,dmax,length(Ns));
sdpt3_times  = zeros(ensemble_size,dmax,length(Ns));
DIA_times    = zeros(ensemble_size,dmax,length(Ns));
sedumi_times = zeros(ensemble_size,dmax,length(Ns));

gdapB_errors    = zeros(ensemble_size,dmax,length(Ns));
mosek_errors    = zeros(ensemble_size,dmax,length(Ns));
sdpt3_errors    = zeros(ensemble_size,dmax,length(Ns));
DIA_errors      = zeros(ensemble_size,dmax,length(Ns));
sedumi_errors   = zeros(ensemble_size,dmax,length(Ns));


for d=2:dmax
    for Nindex=1:length(Npows)
        for i=1:ensemble_size
%             dir = sprintf('./benchmarking_results/d%i',d);
            dir = sprintf('./Ndependence_benchmarking_results/d%i/Npow%i',d,Npows(Nindex));
            clear choi_ground
            load([dir,'/dataset',num2str(i)]);

            clear choi_ml_vec

            load([dir,'/DIA_results',num2str(i)]);
            DIA_times(i,d,Nindex)  = elapsedTime;
            choi_DIA        = reshape(choi_ml_vec,[],d*d);
            DIA_errors(i,d,Nindex) = trace_dist(choi_ground/trace(choi_ground),choi_DIA/trace(choi_DIA));

            clear choi_ml_vec

            load([dir,'/gdapB_results',num2str(i)]);
            gdapB_times(i,d,Nindex) = elapsedTime;
            choi_gdapB        = reshape(choi_ml_vec,[],d*d);
            gdapB_errors(i,d,Nindex) = trace_dist(choi_ground/trace(choi_ground),choi_gdapB/trace(choi_gdapB));

            clear choi_ml_vec

            load([dir,'/mosek_results',num2str(i)]);
            mosek_times(i,d,Nindex) = elapsedTime;    
            choi_mosek        = reshape(choi_ml_vec,[],d*d);
            mosek_errors(i,d,Nindex) = trace_dist(choi_ground/trace(choi_ground),choi_mosek/trace(choi_mosek));

%             clear choi_ml_vec
% 
%             load([dir,'/sdpt3_results',num2str(i)]);
%             sdpt3_times(i,d,Nindex) = elapsedTime;    
%             choi_sdpt3        = reshape(choi_ml_vec,[],d*d);
%             sdpt3_errors(i,d,Nindex) = trace_dist(choi_ground/trace(choi_ground),choi_sdpt3/trace(choi_sdpt3));

%             clear choi_ml_vec
% 
%             load([dir,'/sedumi_results',num2str(i)]);
%             sedumi_times(i,d,Nindex) = elapsedTime;    
%             choi_sdpt3        = reshape(choi_ml_vec,[],d*d);
%             sedumi_errors(i,d,Nindex) = trace_dist(choi_ground/trace(choi_ground),choi_sedumi/trace(choi_sedumi));

        end 
    end
    
    figure; hold on;
    
    errorbar(Ns,squeeze(mean(DIA_times(:,d,:)))',squeeze(std(DIA_times(:,d,:)))','-d','LineWidth',2);
    errorbar(Ns,squeeze(mean(gdapB_times(:,d,:)))',squeeze(std(gdapB_times(:,d,:)))','-*','LineWidth',2);
    errorbar(Ns,squeeze(mean(mosek_times(:,d,:)))',squeeze(std(mosek_times(:,d,:)))','-x','LineWidth',2);
%     errorbar(Ns,squeeze(mean(sdpt3_times(:,d,:)))',squeeze(std(sdpt3_times(:,d,:)))','-x','LineWidth',2);

                    
    xlabel 'N'
    ylabel 'times taken (s)';
    set(gca,'XScale','log')
    set(gca,'YScale','log');
    legend('DIA','gdapB','mosek')
    % legend('gdapB','mosek','sdpt3')
    box on
    grid on
    set(gca,'fontsize',18)
    title(['d = ',num2str(d)])
    saveas(gcf,['./plots/timed',num2str(d),'.png'])
%     saveas(gcf,['./plots/timed',num2str(d),'.eps'],'epsc')
    set(gcf,'paperpositionmode','auto')
    print(gcf,'-depsc2','-loose',['./plots/timed',num2str(d),'.eps'])


    figure; hold on;
    
    errorbar(Ns,squeeze(mean(DIA_errors(:,d,:)))',squeeze(std(DIA_errors(:,d,:)))','-d','LineWidth',2);
    errorbar(Ns,squeeze(mean(gdapB_errors(:,d,:)))',squeeze(std(gdapB_errors(:,d,:)))','-*','LineWidth',2);
    errorbar(Ns,squeeze(mean(mosek_errors(:,d,:)))',squeeze(std(mosek_errors(:,d,:)))','-x','LineWidth',2);
%     errorbar(Ns,squeeze(mean(sdpt3_errors(:,d,:)))',squeeze(std(sdpt3_errors(:,d,:)))','-x','LineWidth',2);
                    
    xlabel 'N'
    ylabel 'error';
  set(gca,'XScale','log');
  set(gca,'YScale','log');
    legend('DIA','gdapB','mosek')
    % legend('gdapB','mosek','sdpt3')
    box on
    grid on
    set(gca,'fontsize',18)
    title(['d = ',num2str(d)])
    saveas(gcf,['./plots/errord',num2str(d),'.png'])
%     saveas(gcf,['./plots/errord',num2str(d),'.eps'],'epsc')
    print(gcf,'-depsc2','-loose',['./plots/errord',num2str(d),'.eps'])
    set(gcf,'paperpositionmode','auto')
end


% figure; hold on;
% scatter(DIA_errors(:,2),DIA_times(:,2),'o','filled','DisplayName','DIA')
% scatter(gdapB_errors(:,2),gdapB_times(:,2),'o','filled','DisplayName','gdapB')
% scatter(mosek_errors(:,2),mosek_times(:,2),'o','filled','DisplayName','mosek')
% scatter(sdpt3_errors(:,2),sdpt3_times(:,2),'o','filled','DisplayName','sdpt3')
% scatter(sedumi_errors(:,2),sedumi_times(:,2),'o','filled','DisplayName','sedumi')
% 
% legend('show','Location','northwest')
% 
% ax = gca;
% ax.ColorOrderIndex = 1;
% 
% scatter(DIA_errors(:,3),DIA_times(:,3),'<','filled')
% scatter(gdapB_errors(:,3),gdapB_times(:,3),'<','filled')
% scatter(mosek_errors(:,3),mosek_times(:,3),'<','filled')
% scatter(sdpt3_errors(:,3),sdpt3_times(:,3),'<','filled')
% 
% ax = gca;
% ax.ColorOrderIndex = 1;
% 
% scatter(DIA_errors(:,4),DIA_times(:,4),'d','filled')
% scatter(gdapB_errors(:,4),gdapB_times(:,4),'d','filled')
% scatter(mosek_errors(:,4),mosek_times(:,4),'d','filled')
% scatter(sdpt3_errors(:,4),sdpt3_times(:,4),'d','filled')
% scatter(sedumi_errors(:,4),sedumi_times(:,4),'d','filled')
% 
% ax = gca;
% ax.ColorOrderIndex = 1;
% 
% scatter(DIA_errors(:,5),DIA_times(:,5),'p','filled')
% scatter(gdapB_errors(:,5),gdapB_times(:,5),'p','filled')
% scatter(mosek_errors(:,5),mosek_times(:,5),'p','filled')
% scatter(sdpt3_errors(:,5),sdpt3_times(:,5),'p','filled')
% scatter(sedumi_errors(:,5),sedumi_times(:,5),'p','filled')
% 
% ax = gca;
% ax.ColorOrderIndex = 1;
% 
% scatter(DIA_errors(:,6),DIA_times(:,6),'s','filled')
% scatter(gdapB_errors(:,6),gdapB_times(:,6),'s','filled')
% scatter(mosek_errors(:,6),mosek_times(:,6),'s','filled')
% scatter(sdpt3_errors(:,6),sdpt3_times(:,6),'s','filled')
% scatter(sedumi_errors(:,6),sedumi_times(:,6),'s','filled')
% 
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% 
% xlabel 'error'
% ylabel 'time (s)'
% 
% box on
% grid on
% 
% 
% 
% set(gca,'fontsize',18)
% saveas(gcf,'scatter.png')
% saveas(gcf,'scatter.eps','epsc')
