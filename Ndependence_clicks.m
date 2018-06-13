% script to investigate algorithm performance as a function of N
% script to generate a number of simulated datasets
% addpath('./QETLAB-0.9')
% addpath('./QETLAB-0.9/helpers')
ensemble_size = 10;
default_ensemble = 'qp';
if exist('ensemble')
else
    ensemble = default_ensemble
end
% figure;
for d=4:4
    
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
    
    
    for l=1:ensemble_size
        fprintf('%d ', l); 
        % generate random ground truth
%         choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
%         choi_ground_vec = reshape(choi_ground,[],1);
%         choi_ground_vec = CPTP_project(choi_ground_vec, MdagM, Mdagb);
%         choi_ground     = reshape(choi_ground_vec,[],d*d);
        
% %         choi_ground     = randomCPTP(d,1); % kraus rank 1, i.e unitary map.
%         choi_ground     = randomCPTP(d,d*d); % kraus rank is full.
%         choi_ground_vec = reshape(choi_ground,[],1);

        switch ensemble
            case 'qp'
                choi_ground     = randomCPTP_quasi_pure(d,0.9);
            case 'fr'
                choi_ground     = randomCPTP(d,d*d); % kraus rank is full.
        end
                
%         partial_trace(choi_ground)
        choi_ground_vec = reshape(choi_ground,[],1);

        p               = real(A*choi_ground_vec);
%         p               = p/sum(p);
%         p               = p*d*d;
        
        for Npow=[1,2,3,4,5,6,7,8,Inf] % above Npow=9 the memory requirements are huge for simulating multinomial noise                            
%         for Npow=[Inf]
            N = 10^Npow
            
            if isinf(N)
%                 p           = reshape(p,[],1);
                n           = p;
%                 n           = n./sum(n);
            else
                pmat           = reshape(p,[],2*d*d); % need an object with n_measurement_outcomes columns
%                 subplot(3,1,1)
%                 bar3(pmat);
%                 hold on;
%                 for m=1:d*d % loop over input states
%                     pmat(m,:) = pmat(m,:)./sum(pmat(m,:));
%                 end
                nmat        = mnrnd(N,pmat);
                nmat        = nmat./sum(nmat,2); % proper normalisation so that sum(nmat,2)=1
                n           = reshape(nmat,[],1);
                
%                 n           = n./sum(n);
%                 sum(n)

        
% 
%                 for i=1:2*d*d:(2*d*d*d*d-d*d) % for each preparation take the click distribution
%                     n(i:i+2*d*d-1) = n(i:i+2*d*d-1)/sum(n(i:i+2*d*d-1)); % normalise to 'frequencies' 
%                     % TODO make sure that this is simply dividing whole
%                     % likelihood by a constant.
%                 end
%                 if Npow ==2
%                     l = 2
%                 else
%                     l = 3
%                 end
%                 subplot(3,1,l)
%                 bar3(reshape(n,[],2*d*d));
            end
%         n               = p; % noiseless scenario

% save A as well, or assume fixed?
            dir = sprintf('./Ndependence_benchmarking_results/d%i/Npow%i',d,Npow);
            save([dir,'/dataset',num2str(l)],'choi_ground','n','p')
            
        end
    end
    
end