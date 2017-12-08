% script to generate a number of simulated datasets
% for varying dimension and using the GGM I/O matrix
addpath('./QETLAB-0.9')
addpath('./QETLAB-0.9/helpers')
ensemble_size = 10;

for d=2:6
    dir = sprintf('./benchmarking_results/d%i',d);
    fprintf(char(10));
    fprintf('%d ', d);
    A = GGM_IO(d);
    
    
    for i=1:ensemble_size
        fprintf('%d ', i); 
        % generate random ground truth
        choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
        choi_ground_vec = reshape(choi_ground,[],1);
        choi_ground_vec = CPTP_project(choi_ground_vec);
        choi_ground     = reshape(choi_ground_vec,[],d*d);

        p               = real(A*choi_ground_vec);
        p               = p/sum(p);
        % n             = mnrnd(1e4,p)';
        % n             = n/sum(n); % activate for multinomial noise
        
        n               = p; % noiseless scenario

% save A as well, or assume fixed?
       
        save([dir,'/dataset',num2str(i)],'choi_ground','n','p')
    end
    
end