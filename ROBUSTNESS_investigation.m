% script to generate a number of simulated datasets

clear;
ensemble_size = 1000;
close all;
d=2;
fig1 = figure;
fig2 = figure;
for rank=1:4
    fprintf('%d: ', d);
    A = PM_minimal(d);  
    
    for i=1:ensemble_size
        fprintf('%d ', i); 
        % generate random ground truth
%         
        choi_ground     = randomCPTP(d,rank);
        choi_ground_vec = reshape(choi_ground,[],1);

        p               = real(A*choi_ground_vec);        
    
        n               = p; % noiseless scenario

        [choi_ml_vecROBUST, ~, ~]  = gdapB(A,n);
        [choi_ml_vecROBUST_IT, ~, ~]  = gdapB_IT(A,n);
        choi_mlROBUST = reshape(choi_ml_vecROBUST,[],d*d);
        choi_mlROBUST_IT = reshape(choi_ml_vecROBUST_IT,[],d*d);

        errorROBUST(i) = norm(choi_mlROBUST-choi_ground);
        errorROBUST_IT(i) = norm(choi_mlROBUST_IT-choi_ground);
        
    end
    figure(fig1); plot(errorROBUST,'LineWidth',2); hold on;
    figure(fig2); plot(errorROBUST_IT,'LineWidth',2); hold on;
end

figure(fig1)
pbaspect([4,1,1])
ylim([1e-7,1])
title('Conditioned PGDB')
legend('rank 1','rank2','rank3','rank4')
set(gca,'YScale','log');
xlabel('run')
ylabel('error')
grid on
box on
set(gca,'fontsize',18)
print(gcf,'-depsc2','-loose',['./plots/ROBUSTNESSd',num2str(d),'.eps'])

figure(fig2)
pbaspect([4,1,1])
ylim([1e-7,1])
title('identity trick')
legend('rank 1','rank2','rank3','rank4')
set(gca,'YScale','log');
xlabel('run')
ylabel('error')
grid on
box on
set(gca,'fontsize',18)
print(gcf,'-depsc2','-loose',['./plots/ROBUSTNESSd_IT',num2str(d),'.eps'])
