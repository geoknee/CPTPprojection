% script to generate a number of simulated datasets
addpath('./QETLAB-0.9')
addpath('./QETLAB-0.9/helpers')
ensemble_size = 30;
default_ensemble = 'qp';
if exist('ensemble')
else
    ensemble = default_ensemble
end

for d=2:5
    dir = sprintf('./benchmarking_results/d%i',d);
    fprintf(char(10));
    fprintf('%d ', d);
    A = PM_minimal(d);
%     A = GGMall_IO(d);
    
    % precompute matrices for TP_project
    M = zeros([d*d,d*d*d*d]);
    for i=1:d
        e = zeros(1,d);
        e(i)  = 1;
        B = kron(eye(d),e); 
        M = M + kron(B,B);
    end
    MdagM = sparse(M'*M);
    b = sparse(reshape(eye(d),[],1));
    Mdagb = sparse(M'*b);
    
    
    for i=1:ensemble_size
        fprintf('%d ', i); 
        % generate random ground truth
%         choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
%         choi_ground_vec = reshape(choi_ground,[],1);
%         choi_ground_vec = CPTP_project(choi_ground_vec, MdagM, Mdagb);
%         choi_ground     = reshape(choi_ground_vec,[],d*d);
        
%         choi_ground     = randomCPTP(d,1); % kraus rank 1, i.e unitary map.
%         choi_ground     = randomCPTP(d,d*d); % kraus rank is full.
        
        switch ensemble
            case 'qp'
                choi_ground     = randomCPTP_quasi_pure(d,0.9);
            case 'fr'
                choi_ground     = randomCPTP(d,d*d); % kraus rank is full.
        end
        
        choi_ground_vec = reshape(choi_ground,[],1);
        
        p               = real(A*choi_ground_vec);
%         p               = p/sum(p);
        % n             = mnrnd(1e4,p)';
        % n             = n/sum(n); % activate for multinomial noise
        
        n               = p; % noiseless scenario
%         n               = n./sum(n); % noiseless scenario

% save A as well, or assume fixed?
       
        save([dir,'/dataset',num2str(i)],'choi_ground','n','p')
    end
    
end