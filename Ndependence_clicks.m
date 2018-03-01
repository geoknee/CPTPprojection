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
    
    
    for l=1:ensemble_size
        fprintf('%d ', l); 
        % generate random ground truth
        choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
        choi_ground_vec = reshape(choi_ground,[],1);
        choi_ground_vec = CPTP_project(choi_ground_vec, MdagM, Mdagb);
        choi_ground     = reshape(choi_ground_vec,[],d*d);

        p               = real(A*choi_ground_vec);
        p               = reshape(p,[],d*d);
%         p               = p/sum(p);
        
        for Npow=[0,1,2,3,4,5,6,7]
            N = 10^Npow;
            n           = reshape(mnrnd(N,p')',[],1);
        
            
            for i=1:2*d*d:(2*d*d*d*d-d*d) % for each preparation take the click distribution
                n(i:i+2*d*d-1) = n(i:i+2*d*d-1)/sum(n(i:i+2*d*d-1)); % normalise to 'frequencies'
            end
%         n               = p; % noiseless scenario

% save A as well, or assume fixed?
            dir = sprintf('./Ndependence_benchmarking_results/d%i/Npow%i',d,Npow);
            save([dir,'/dataset',num2str(l)],'choi_ground','n','p')
            
        end
    end
    
end