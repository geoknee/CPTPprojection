% script to investigate algorithm performance as a function of N
% script to generate a number of simulated datasets
% addpath('./QETLAB-0.9')
% addpath('./QETLAB-0.9/helpers')
ensemble_size = 10;

for d=2:5
    
    fprintf(newline);
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
        choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
        choi_ground_vec = reshape(choi_ground,[],1);
        choi_ground_vec = CPTP_project(choi_ground_vec, MdagM, Mdagb);
        choi_ground     = reshape(choi_ground_vec,[],d*d);

        p               = real(A*choi_ground_vec);
        p               = p/sum(p);
        
        for N=[2,4,8,16,32,64,128,256,512,1024,2048,4096]
            n             = mnrnd(N,p)';
            n             = n/sum(n); %  multinomial noise
        
%         n               = p; % noiseless scenario

% save A as well, or assume fixed?
            dir = sprintf('./Ndependence_benchmarking_results/d%i/N%i',d,N);
            save([dir,'/dataset',num2str(i)],'choi_ground','n','p')
            
        end
    end
    
end