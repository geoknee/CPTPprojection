% use a varitey of algorithms to find the maximum likelihood process 
% from various datasets 
% parpool(3)
clear all
% addpath('./QETLAB-0.9')
% addpath('./QETLAB-0.9/helpers')
ensemble_size = 1;

for d=4:4   
%     for method={'mosek','gdapB','DIA'}
for method={'gdapB'}
% for method = {'DIA'}
        fprintf(char(10));
        fprintf(method{1});

        fprintf(char(10));
        fprintf('%d :', d);

        A = PM_minimal(d);
%          A = GGMall_IO(d);
        for Npow=[1,3,5,7,Inf]
            N = 10^Npow;
            dir = sprintf('./Ndependence_benchmarking_results/d%i/Npow%i',d,Npow);
            fprintf(char(10));
            fprintf('[%d] ', N);
            for i = 1:ensemble_size % total number of simulated datasets
                fprintf('%d ', i); 
                load([dir,'/dataset',num2str(i)]);

                    switch char(method)
                        case'gdapB'
                            tic;
                            [choi_ml_vec, solution, costs] = gdapB(A,n);
                            elapsedTime = toc;
    %                         semilogy(costs)
    %                         hold on
                        case 'DIA'
                            tic;
                            [choi_ml_vec, solution, costs] = DIA(A,n);
                            elapsedTime = toc
                            semilogy(costs)
                            hold on
                        case 'mosek'
                            cvx_solver mosek
                            tic;
                            [choi_ml_vec] = cvx_wrapper(A,n);
                            elapsedTime = toc;
                            solution=[]; costs = []; % cannot currently extract these
                        case 'sdpt3'
                            cvx_solver sdpt3
                            tic;
                            [choi_ml_vec] = cvx_wrapper(A,n);
                            elapsedTime = toc;
                            solution=[]; costs = []; % cannot currently extract these
                        case 'sedumi'
                            cvx_solver sedumi
                            tic;
                            [choi_ml_vec] = cvx_wrapper(A,n);
                            elapsedTime = toc;
                            solution=[]; costs = []; % cannot currently extract these
                    end
                    save([dir,'/',char(method),'_results',num2str(i)],'elapsedTime','choi_ml_vec','costs','solution')
            end
        end
    end
end


