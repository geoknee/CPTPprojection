% read in ensemble_size running times and precisions for each method, for each d.
%% check global variable set
if exist('ensemble')
   fprintf(['ensemble = ',ensemble])
else
    error('you must set the ensemble variable to either qp (quasi pure) or fr (full rank)')
end

if exist('drange')
    fprintf(['drange = ',drange])
else
    error('you must set the drange variable')
end

if exist('LIswitch')
    fprintf(['LIswitch = ',LIswitch])
else
    error('you must set the LIsiwtch variable (if 0 Linear Inversion is run, if 1 it is not)')
end

if exist('ensemble_size')
    fprintf(['ensemble_size = ',ensemble_size])
else
    error('you must set the ensemble_size variable')
end

if exist('Npows')
    fprintf(['Npows = ',ensemble_size])
else
    error('you must set the Npows variable')
end
%%
Ns = 10.^Npows;
Ns(end)=10^12;
dmax = drange(end);
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


for d=drange
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

            clear choi_ml_vec
            
            load([dir,'/LinInversion_results',num2str(i)]);
            LI_times(i,d,Nindex) = elapsedTime;    
            choi_LI        = reshape(choi_ml_vec,[],d*d);
            LI_errors(i,d,Nindex) = trace_dist(choi_ground/trace(choi_ground),choi_LI/trace(choi_LI));

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
    %%
    figure; ax1 = subplot(2,2,1); hold on;
    
    errorbar(Ns(1:end-1),squeeze(nanmean(DIA_times(:,d,1:end-1)))',squeeze(nanstd(DIA_times(:,d,1:end-1)))','-d','LineWidth',2);
    errorbar(Ns(1:end-1),squeeze(nanmean(gdapB_times(:,d,1:end-1)))',squeeze(nanstd(gdapB_times(:,d,1:end-1)))','-*','LineWidth',2);
    errorbar(Ns(1:end-1),squeeze(nanmean(mosek_times(:,d,1:end-1)))',squeeze(nanstd(mosek_times(:,d,1:end-1)))','-x','LineWidth',2);
    
    if LIswitch
        errorbar(Ns(1:end-1),squeeze(nanmean(LI_times(:,d,1:end-1)))',squeeze(nanstd(LI_times(:,d,1:end-1)))','-x','LineWidth',2);
    end

    %     errorbar(Ns,squeeze(nanmean(sdpt3_times(:,d,:)))',squeeze(nanstd(sdpt3_times(:,d,:)))','-x','LineWidth',2);

                    
%     xlabel 'N'
    ylabel 'times taken (s)';
    set(gca,'XScale','log')
%     set(gca,'YScale','log');
    set(gca,'XTick',Ns(1:end-1))
%     set(gca,'YTick',[1e-3,1e-2,1e-1,1e0,1e1,1e2])
    box on; grid on;
    set(gca,'XTickLabel',[])
%     legend('DIA','gdapB','mosek','Location','NorthWest')
    % legend('gdapB','mosek','sdpt3')
    box on
    grid off
%     set(gca,'fontsize',18)
    
    ax2 = subplot(2,2,2); hold on; pbaspect([1 4 1]);
    
    errorbar(Ns(end),squeeze(nanmean(DIA_times(:,d,end)))',squeeze(nanstd(DIA_times(:,d,end)))','-d','LineWidth',2);
    errorbar(Ns(end),squeeze(nanmean(gdapB_times(:,d,end)))',squeeze(nanstd(gdapB_times(:,d,end)))','-*','LineWidth',2);
    errorbar(Ns(end),squeeze(nanmean(mosek_times(:,d,end)))',squeeze(nanstd(mosek_times(:,d,end)))','-x','LineWidth',2);
    if LIswitch
        errorbar(Ns(end),squeeze(nanmean(LI_times(:,d,end)))',squeeze(nanstd(LI_times(:,d,end)))','-x','LineWidth',2);
    end
    
    set(gca,'XScale','log');
%     set(gca,'YScale','log');
    box on; grid off;
    ax = gca;
    set(gca,'XTick',[])
    ax.YAxisLocation = 'right';
%     title(['d = ',num2str(d)])
%     saveas(gcf,['./plots/timed',num2str(d),'.png'])
% %     saveas(gcf,['./plots/timed',num2str(d),'.eps'],'epsc')
%     set(gcf,'paperpositionmode','auto')
%     print(gcf,'-depsc2','-loose',['./plots/timed',num2str(d),'.eps'])


%     figure; 
    ax3 = subplot(2,2,3); hold on;  
    
    errorbar(Ns(1:end-1),squeeze(nanmean(DIA_errors(:,d,1:end-1)))',squeeze(nanstd(DIA_errors(:,d,1:end-1)))','-d','LineWidth',2);
    errorbar(Ns(1:end-1),squeeze(nanmean(gdapB_errors(:,d,1:end-1)))',squeeze(nanstd(gdapB_errors(:,d,1:end-1)))','-*','LineWidth',2);
    errorbar(Ns(1:end-1),squeeze(nanmean(mosek_errors(:,d,1:end-1)))',squeeze(nanstd(mosek_errors(:,d,1:end-1)))','-x','LineWidth',2);
    if LIswitch
        errorbar(Ns(1:end-1),squeeze(nanmean(LI_errors(:,d,1:end-1)))',squeeze(nanstd(LI_errors(:,d,1:end-1)))','-x','LineWidth',2);
    end
    %     errorbar(Ns,squeeze(nanmean(sdpt3_errors(:,d,:)))',squeeze(nanstd(sdpt3_errors(:,d,:)))','-x','LineWidth',2);


    xlabel 'N'
    ylabel 'error';
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    set(gca,'XTick',[1e1,1e3,1e5,1e7,1e9])
    set(gca,'YTick',[1e-5,1e-4,1e-3,1e-2,1e-1,1e0])
    box on; grid on;
    if LIswitch
        legend('DIA','pgdB','mosek','LIFP','Location','southwest')
    else
        legend('DIA','pgdB','mosek','Location','southwest')
    end
    
    ax4 = subplot(2,2,4); hold on;pbaspect([1 4 1]);
    errorbar(Ns(end),squeeze(nanmean(DIA_errors(:,d,end)))',squeeze(nanstd(DIA_errors(:,d,end)))','-d','LineWidth',2);
    errorbar(Ns(end),squeeze(nanmean(gdapB_errors(:,d,end)))',squeeze(nanstd(gdapB_errors(:,d,end)))','-*','LineWidth',2);
    errorbar(Ns(end),squeeze(nanmean(mosek_errors(:,d,end)))',squeeze(nanstd(mosek_errors(:,d,end)))','-x','LineWidth',2);
    if LIswitch    
        errorbar(Ns(end),squeeze(nanmean(LI_errors(:,d,end)))',squeeze(nanstd(LI_errors(:,d,end)))','-x','LineWidth',2);
    end
    
    xlabel('$\infty$','Interpreter','latex')
    set(gca,'YTick',[1e-5,1e-4,1e-3,1e-2,1e-1,1e0])
    set(gca,'XTick',[])
%     ylabel 'error';
    ax=gca;
    ax.YAxisLocation = 'right';
    set(gca,'XScale','log');
    set(gca,'YScale','log');
%     legend('DIA','gdapB','mosek')
    % legend('gdapB','mosek','sdpt3')
    
    
    box on
    grid on
%     set(gca,'fontsize',18)
%     title(['d = ',num2str(d)])
   
    linkaxes([ax1,ax2],'y')
%     linkaxes([ax3,ax4],'y')
%     linkaxes([ax1,ax3],'x')
%     linkaxes([ax2,ax4],'x')

    saveas(gcf,['./plots/',ensemble,'timeerrord',num2str(d),'.png'])
    saveas(gcf,['./plots/',ensemble,'ALTtimeerrord',num2str(d),'.eps'],'epsc')
    set(gcf,'paperpositionmode','auto')
    print(gcf,'-depsc2','-loose',['./plots/',ensemble,'timeerrord',num2str(d),'.eps'])
    savefig(['./plots/',ensemble,'timeerrord',num2str(d)])
end

