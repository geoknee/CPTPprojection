% script to investigate algorithm performance as a function of N
% script to generate a number of simulated datasets
% addpath('./QETLAB-0.9')
% addpath('./QETLAB-0.9/helpers')
ensemble_size = 10;

for d=4:4
    
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
%         choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
%         choi_ground_vec = reshape(choi_ground,[],1);
%         choi_ground_vec = CPTP_project(choi_ground_vec, MdagM, Mdagb);
%         choi_ground     = reshape(choi_ground_vec,[],d*d);
        
%         choi_ground     = randomCPTP(d,1); % kraus rank 1, i.e unitary map.
        choi_ground     = randomCPTP(d,d*d); % kraus rank is full.
        choi_ground_vec = reshape(choi_ground,[],1);

        p               = real(A*choi_ground_vec);
        
%         p               = p/sum(p);
        
        for Npow=[1,2,3,4,5,6,7,8,9,Inf] % above Npow=9 the memory requirements are huge for simulating multinomial noise
                                    
            N = 10^Npow;
            
            if isinf(N)
                p           = reshape(p,[],1);
                n           = p;
            else
                p           = reshape(p,[],d*d);
                n           = reshape(mnrnd(N,p')',[],1);
        

                for i=1:2*d*d:(2*d*d*d*d-d*d) % for each preparation take the click distribution
                    n(i:i+2*d*d-1) = n(i:i+2*d*d-1)/sum(n(i:i+2*d*d-1)); % normalise to 'frequencies'
                end
            end
%         n               = p; % noiseless scenario

% save A as well, or assume fixed?
            dir = sprintf('./Ndependence_benchmarking_results/d%i/Npow%i',d,Npow);
            save([dir,'/dataset',num2str(l)],'choi_ground','n','p')
            
        end
    end
    
end