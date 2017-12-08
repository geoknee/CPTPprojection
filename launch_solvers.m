% run 
% parpool(3)
addpath('./QETLAB-0.9')
addpath('./QETLAB-0.9/helpers')
ensemble_size = 10;

cvx_solver mosek
for method={'mosek','gdapB','DIA'}
% for method={'gdapB'}
    for d=2:6
        fprintf(char(10));
        fprintf('%d ', d);

        dir = sprintf('./benchmarking_results/d%i',d);

        A = GGM_IO(d);
        for i = 1:ensemble_size % total number of simulated datasets
            fprintf('%d ', i); 
            load([dir,'/dataset',num2str(i)]);

                switch char(method)
                    case'gdapB'
                        tic;
                        [choi_ml_vec,solution, costs] = gdapB(A,n);
                        elapsedTime = toc;
                    case 'DIA'
                        tic;
                        [choi_ml_vec, solution, costs] = DIA(A,n);
                        elapsedTime = toc;
                    case 'mosek'  
                        tic;
                        [choi_ml_vec] = mosek(A,n);
                        elapsedTime = toc;
                        solution=[]; costs = []; % cannot currently extract these
                end
                save([dir,'/',char(method),'_results',num2str(i)],'elapsedTime','choi_ml_vec','costs','solution')
        end
    end
end



