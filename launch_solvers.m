% run 
% parpool(3)

cvx_solver mosek
for d=2:6
    fprintf(char(10));
    fprintf('%d ', d);
    
    dir = sprintf('./benchmarking_results/d%i',d);

    A = GGM_IO(d);

    for i = 1:100 % total number of simulated datasets
        fprintf('%d ', i); 
        load([dir,'/dataset',num2str(i)]);
        
        for method={'gdapB','DIA','mosek'}
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



