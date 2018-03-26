% script to generate a number of simulated datasets
% addpath('./QETLAB-0.9')
% addpath('./QETLAB-0.9/helpers')
clear;
ensemble_size = 100;
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
%         choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
%         choi_ground_vec = reshape(choi_ground,[],1);
%         choi_ground_vec = CPTNI_project(choi_ground_vec);
%         choi_ground     = reshape(choi_ground_vec,[],d*d);
%         
        choi_ground     = randomCPTP(d,rank); % kraus rank 1, i.e unitary map.
%         choi_ground     = randomCPTP(d,d*d); % kraus rank is full.
        choi_ground_vec = reshape(choi_ground,[],1);
        
%         
%         if min(eig(choi_ground))<-1e-16
%             sprintf('choi matrix not PSD')
%             eig(choi_ground)
%         end

        p               = real(A*choi_ground_vec);        
    
        n               = p; % noiseless scenario

%         N           = 1e4;
%         p           = reshape(p,[],d*d);
%         for j = 1:d*d
%             p(2*d*d+1,j)  = 1-sum(p(:,j)); % expand to a TP process in d+1
%         end
%         n           = mnrnd(N,p')';
%         n           = n(1:2*d*d,:); % discard binned state
%         n           =reshape(n,[],1);

%       
%         [choi_ml_vec, ~, ~] = gdapB(A,n);
%         choi_ml = reshape(choi_ml_vec,[],d*d);
        [choi_ml_vecROBUST, ~, ~]  = gdapB(A,n);
        [choi_ml_vecROBUST_IT, ~, ~]  = gdapB_IT(A,n);
        choi_mlROBUST = reshape(choi_ml_vecROBUST,[],d*d);
        choi_mlROBUST_IT = reshape(choi_ml_vecROBUST_IT,[],d*d);

% 
%         errorTP(i) = trace_dist(choi_mlTP/trace(choi_mlTP),choi_ground/trace(choi_ground));
%         errorTNI(i) = trace_dist(choi_mlTNI/trace(choi_mlTNI),choi_ground/trace(choi_ground)); % is this a good measure for TNI processes?
        
%         error(i) = norm(choi_ml-choi_ground);
        errorROBUST(i) = norm(choi_mlROBUST-choi_ground);
        errorROBUST_IT(i) = norm(choi_mlROBUST_IT-choi_ground);
        
% % save A as well, or assume fixed?

%         save([dir,'/dataset',num2str(i)],'choi_ground','n','p')
    
    end
    figure(fig1); plot(errorROBUST,'LineWidth',2); hold on;
    figure(fig2); plot(errorROBUST_IT,'LineWidth',2); hold on;
end

figure(fig1)
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
title('identity trick')
legend('rank 1','rank2','rank3','rank4')
set(gca,'YScale','log');
xlabel('run')
ylabel('error')
grid on
box on
set(gca,'fontsize',18)
print(gcf,'-depsc2','-loose',['./plots/ROBUSTNESSd_IT',num2str(d),'.eps'])
