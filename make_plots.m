% read in ensemble_size running times and precisions for each method, for each d.
clear;close all;
dmax = 4;
ensemble_size =20;
gdapB_times = zeros(ensemble_size,dmax);
mosek_times = zeros(ensemble_size,dmax);
sdpt3_times = zeros(ensemble_size,dmax);
DIA_times   = zeros(ensemble_size,dmax);

gdapB_errors = zeros(ensemble_size,dmax);
mosek_errors = zeros(ensemble_size,dmax);
sdpt3_errors = zeros(ensemble_size,dmax);
DIA_errors   = zeros(ensemble_size,dmax);


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
        
        
    end
end

figure; hold on;

errorbar(2:dmax,mean(DIA_times(:,2:end)),std(DIA_times(:,2:end)),'-rd');
errorbar(2:dmax,mean(gdapB_times(:,2:end)),std(gdapB_times(:,2:end)),'-b*');
errorbar(2:dmax,mean(mosek_times(:,2:end)),std(mosek_times(:,2:end)),'-gx');
errorbar(2:dmax,mean(sdpt3_times(:,2:end)),std(sdpt3_times(:,2:end)),'-kx');
xlim([1.8,dmax+0.2])
xlabel 'Hilbert space dimension'
ylabel 'times taken (s)';
set(gca,'YScale','log');
legend('DIA','gdapB','mosek','sdpt3')
saveas(gcf,'time.png')


figure; hold on;

errorbar(2:dmax,mean(DIA_errors(:,2:end)),std(DIA_errors(:,2:end)),'-rd');
errorbar(2:dmax,mean(gdapB_errors(:,2:end)),std(gdapB_errors(:,2:end)),'-b*');
errorbar(2:dmax,mean(mosek_errors(:,2:end)),std(mosek_errors(:,2:end)),'-gx');
errorbar(2:dmax,mean(sdpt3_errors(:,2:end)),std(sdpt3_errors(:,2:end)),'-kx');
xlim([1.8,dmax+0.2])
ylim([0,1])
xlabel 'Hilbert space dimension'
ylabel 'J distance';
set(gca,'YScale','log')
legend('DIA','gdapB','mosek','sdpt3')
saveas(gcf,'errors.png')


