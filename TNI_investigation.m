% script to generate a number of simulated datasets
% addpath('./QETLAB-0.9')
% addpath('./QETLAB-0.9/helpers')
clear;
ensemble_size = 10;


for d=2:2
    dir = sprintf('./TNIbenchmarking_results/d%i',d);
    fprintf(newline);
    fprintf('%d: ', d);
    A = PM_minimal(d);  
    
    for i=1:ensemble_size
        fprintf('%d ', i); 
        % generate random ground truth
        choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
        choi_ground_vec = reshape(choi_ground,[],1);
        choi_ground_vec = CPTNI_project(choi_ground_vec);
        choi_ground     = reshape(choi_ground_vec,[],d*d);
        
%         choi_ground     = randomCPTP(d,1); % kraus rank 1, i.e unitary map.
%         choi_ground     = randomCPTP(d,d*d); % kraus rank is full.
%         choi_ground_vec = reshape(choi_ground,[],1);
        
        
        if min(eig(choi_ground))<-1e-16
            sprintf('choi matrix not PSD')
            eig(choi_ground)
        end

        p               = real(A*choi_ground_vec);        
    
        dn               = p; % noiseless scenario

%         N           = 1e4;
%         p           = reshape(p,[],d*d);
%         for j = 1:d*d
%             p(2*d*d+1,j)  = 1-sum(p(:,j)); % expand to a TP process in d+1
%         end
%         n           = mnrnd(N,p')';
%         n           = n(1:2*d*d,:); % discard binned state
%         n           =reshape(n,[],1);

%       
        [choi_ml_vecTNI, ~, ~] = TNIgdapB(A,n);
        choi_mlTNI = reshape(choi_ml_vecTNI,[],d*d);

        [choi_ml_vecTP, ~, ~]  = gdapB(A,n);
        choi_mlTP = reshape(choi_ml_vecTP,[],d*d);


% 
%         errorTP(i) = trace_dist(choi_mlTP/trace(choi_mlTP),choi_ground/trace(choi_ground));
%         errorTNI(i) = trace_dist(choi_mlTNI/trace(choi_mlTNI),choi_ground/trace(choi_ground)); % is this a good measure for TNI processes?
        
        errorTP(i) = norm(choi_mlTP-choi_ground);
        errorTNI(i) = norm(choi_mlTNI-choi_ground);
        
% % save A as well, or assume fixed?

%         save([dir,'/dataset',num2str(i)],'choi_ground','n','p')
    end
            figure;plot(errorTP);hold on; plot(errorTNI); legend('TP assumed','TNI assumed')
end