% script to generate a number of simulated datasets
% addpath('./QETLAB-0.9')
% addpath('./QETLAB-0.9/helpers')

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
for d=drange
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
    b = sparse(reshape(speye(d),[],1));
    Mdagb = sparse(M'*b);
    
    
    for i=1:ensemble_size
        fprintf('%d ', i); 
        
        % generate random ground truth
               
        switch ensemble
            case 'qp'
                choi_ground     = randomCPTP_quasi_pure(d,0.9);
%                 choi_ground     = randomCPTP(d,1); % kraus rank 1, i.e unitary map.

            case 'fr'
                choi_ground     = randomCPTP(d,d*d); % kraus rank is full.
        end
        
        choi_ground_vec = reshape(choi_ground,[],1);
        
        p               = real(A*choi_ground_vec);
       
        n               = p; % noiseless scenario
       
        save([dir,'/dataset',num2str(i)],'choi_ground','n','p')
    end
    
end