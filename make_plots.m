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


%%
dmax = drange(end);

gdapB_times  = zeros(ensemble_size,dmax);
gdapM_times  = zeros(ensemble_size,dmax);
mosek_times  = zeros(ensemble_size,dmax);
sdpt3_times  = zeros(ensemble_size,dmax);
DIA_times    = zeros(ensemble_size,dmax);
% sedumi_times = zeros(ensemble_size,dmax); % this fails
li_times     = zeros(ensemble_size,dmax);

gdapB_errors    = zeros(ensemble_size,dmax);
gdapM_errors    = zeros(ensemble_size,dmax);
mosek_errors    = zeros(ensemble_size,dmax);
sdpt3_errors    = zeros(ensemble_size,dmax);
DIA_errors      = zeros(ensemble_size,dmax);
li_errors       = zeros(ensemble_size,dmax);
% sedumi_errors   = zeros(ensemble_size,dmax);


for d=drange
    d
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
        
        
%         load([dir,'/gdapM_results',num2str(i)]);
%         gdapM_times(i,d) = elapsedTime;
%         choi_gdapM        = reshape(choi_ml_vec,[],d*d);
%         gdapM_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_gdapM/trace(choi_gdapM));
%         
%         clear choi_ml_vec
        
        
        load([dir,'/mosek_results',num2str(i)]);
        mosek_times(i,d) = elapsedTime;    
        choi_mosek        = reshape(choi_ml_vec,[],d*d);
        mosek_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_mosek/trace(choi_mosek));
        

%         
%         load([dir,'/sdpt3_results',num2str(i)]);
%         sdpt3_times(i,d) = elapsedTime;    
%         choi_sdpt3        = reshape(choi_ml_vec,[],d*d);
%         sdpt3_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_sdpt3/trace(choi_sdpt3));
        
%         clear choi_ml_vec
%         
%         load([dir,'/sedumi_results',num2str(i)]);
%         sedumi_times(i,d) = elapsedTime;    
%         choi_sedumi        = reshape(choi_ml_vec,[],d*d);
%         sedumi_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_sedumi/trace(choi_sedumi));

        if LIswitch
            clear choi_ml_vec
            load([dir,'/LinInversion_results',num2str(i)]);
            li_times(i,d)  = elapsedTime;    
            choi_li        = reshape(choi_ml_vec,[],d*d);
            li_errors(i,d) = trace_dist(choi_ground/trace(choi_ground),choi_li/trace(choi_li));
        end
        
    end
end

figure; hold on;

errorbar(2:dmax,nanmean(DIA_times(:,2:end)),nanstd(DIA_times(:,2:end)),'-d','LineWidth',2);
errorbar(2:dmax,nanmean(gdapB_times(:,2:end)),nanstd(gdapB_times(:,2:end)),'-*','LineWidth',2);
% errorbar(2:dmax,nanmean(gdapM_times(:,2:end)),nanstd(gdapB_times(:,2:end)),'-*','LineWidth',2);
errorbar(2:dmax,nanmean(mosek_times(:,2:end)),nanstd(mosek_times(:,2:end)),'-x','LineWidth',2);
% errorbar(2:dmax,nanmean(sdpt3_times(:,2:end)),nanstd(sdpt3_times(:,2:end)),'-x','LineWidth',2);
% errorbar(2:dmax,nanmean(sedumi_times(:,2:end)),nanstd(sedumi_times(:,2:end)),'-s','LineWidth',2);
if LIswitch
    errorbar(2:dmax,nanmean(li_times(:,2:end)),nanstd(li_times(:,2:end)),'-x','LineWidth',2);
end
xlim([1.8,dmax+0.2])
xlabel 'Hilbert space dimension'
ylabel 'times taken (s)';
set(gca,'YScale','log');
legend('DIA','gdapB','mosek')
% legend('gdapB','mosek','sdpt3')
box on
grid on
set(gca,'fontsize',18)
saveas(gcf,'./plots/time.png')
saveas(gcf,'./plots/time.eps','epsc')


figure; hold on;

errorbar(2:dmax,nanmean(DIA_errors(:,2:end)),nanstd(DIA_errors(:,2:end)),'-d','LineWidth',2);
errorbar(2:dmax,nanmean(gdapB_errors(:,2:end)),nanstd(gdapB_errors(:,2:end)),'-*','LineWidth',2);
% errorbar(2:dmax,nanmean(gdapM_errors(:,2:end)),nanstd(gdapB_errors(:,2:end)),'-*','LineWidth',2);
errorbar(2:dmax,nanmean(mosek_errors(:,2:end)),nanstd(mosek_errors(:,2:end)),'-x','LineWidth',2);
% errorbar(2:dmax,nanmean(sdpt3_errors(:,2:end)),nanstd(sdpt3_errors(:,2:end)),'-x','LineWidth',2);
% errorbar(2:dmax,nanmean(sedumi_errors(:,2:end)),nanstd(sedumi_errors(:,2:end)),'-s','LineWidth',2);
if LIswitch
    errorbar(2:dmax,nanmean(li_errors(:,2:end)),nanstd(li_errors(:,2:end)),'-x','LineWidth',2);
end
xlim([1.8,dmax+0.2])
ylim([0,1])
xlabel 'Hilbert space dimension'
ylabel 'J distance';
set(gca,'YScale','log')
legend('DIA','gdapB','mosek')
% legend('gdapB','mosek','sdpt3')
box on
grid on
set(gca,'fontsize',18)
saveas(gcf,'./plots/errors.png')
saveas(gcf,'./plots/errors.eps','epsc')





%%
figure('Position',[1 0 400 250]); hold on;

for d = 2:dmax
    ax = gca;
    ax.ColorOrderIndex = 1;
    ax.YAxisLocation = 'right';
    
    DIAx(d)    = nanmean(DIA_errors(:,d));
    DIAxbar(d) = nanstd(DIA_errors(:,d));
    DIAy(d)    = nanmean(DIA_times(:,d));
    DIAybar(d) = nanstd(DIA_times(:,d));
    
    errorbar(DIAx(d),DIAy(d),DIAybar(d),DIAybar(d),DIAxbar(d),DIAxbar(d),'LineWidth',2)
       
    gdapBx(d)    = nanmean(gdapB_errors(:,d));
    gdapBxbar(d) = nanstd(gdapB_errors(:,d));
    gdapBy(d)    = nanmean(gdapB_times(:,d));
    gdapBybar(d) = nanstd(gdapB_times(:,d));
    
    errorbar(gdapBx(d),gdapBy(d),gdapBybar(d),gdapBybar(d),gdapBxbar(d),gdapBxbar(d),'LineWidth',2)
    
    mosekx(d)    = nanmean(mosek_errors(:,d));
    mosekxbar(d) = nanstd(mosek_errors(:,d));
    moseky(d)    = nanmean(mosek_times(:,d));
    mosekybar(d) = nanstd(mosek_times(:,d));
    
    errorbar(mosekx(d),moseky(d),mosekybar(d),mosekybar(d),mosekxbar(d),mosekxbar(d),'LineWidth',2)
    
    if LIswitch
        lix(d)    = nanmean(li_errors(:,d));
        lixbar(d) = nanstd(li_errors(:,d));
        liy(d)    = nanmean(li_times(:,d));
        liybar(d) = nanstd(li_times(:,d));
        errorbar(lix(d),liy(d),liybar(d),liybar(d),lixbar(d),lixbar(d),'LineWidth',2)
    end

    
    
    %triangles
%     plot([DIAx(d),gdapBx(d)],[DIAy(d),gdapBy(d)],'k')
%     plot([DIAx,mosekx],[DIAy,moseky],'k')
%     plot([mosekx(d),gdapBx(d)],[moseky(d),gdapBy(d)],'k')

end
ax = gca;
ax.ColorOrderIndex = 1;
   
plot(DIAx,DIAy,'LineWidth',2)
plot(gdapBx,gdapBy,'LineWidth',2)
plot(mosekx,moseky,'LineWidth',2)
if LIswitch
    plot(lix,liy,'LineWidth',2)
end


set(gca,'xscale','log')
set(gca,'yscale','log')

xlabel 'error'
ylabel 'time (s)'

set(gca,'XTick',([1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0]))
% xlim([5e-5,5e-2])
ylim([1e-2,1e2])
box on
grid on
legend('DIA','pgdB','mosek','Location','northwest')

set(gca,'fontsize',12)
saveas(gcf,['./plots/',ensemble,'triangles.png'])
saveas(gcf,['./plots/',ensemble,'triangles.eps'],'epsc')
savefig(['./plots/',ensemble,'triangles'])
