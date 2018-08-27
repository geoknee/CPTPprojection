% use a varitey of algorithms to perform process tomography
% from various datasets 
% parpool(3)
% clear all
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
    error('you must set the LIsiwtch variable (if 1 Linear Inversion is run, if 0 it is not)')
end

if exist('ensemble_size')
    fprintf(['ensemble_size = ',ensemble_size])
else
    error('you must set the ensemble_size variable')
end


%%
if LIswitch
    list_of_methods = {'mosek','gdapB','DIA','LinInversion'};
else
    list_of_methods = {'mosek','gdapB','DIA'};
end

for d=drange
    for method=list_of_methods

        fprintf(char(10));
        fprintf(method{1});

        fprintf(char(10));
        fprintf('%d :', d);

        dir = sprintf('./benchmarking_results/d%i',d);

         A = PM_minimal(d);
%          A = GGMall_IO(d);
      
        for i = 1:ensemble_size % total number of simulated datasets
            fprintf('%d ', i); 
            load([dir,'/dataset',num2str(i)]);
%             n = n./sum(n);

                switch char(method)
                    case'gdapM'
                        tic;
                        [choi_ml_vec, solution, costs] = gdapM(A,n);
                        elapsedTime = toc;
                    case'gdapB'
                        tic;
                        [choi_ml_vec, solution, costs] = gdapB(A,n);
                        elapsedTime = toc;
%                         semilogy(costs)
%                         hold on
                    case 'DIA'
                        tic;
                        [choi_ml_vec, solution, costs] = DIA(A,n);
                        elapsedTime = toc;
%                         semilogy(costs)
%                         hold on
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
                    case 'LinInversion'
                        tic;
                        [choi_ml_vec] = LinInversion(A,n);
                        elapsedTime = toc;
                        solution = []; costs = [];
                end
                save([dir,'/',char(method),'_results',num2str(i)],'elapsedTime','choi_ml_vec','costs','solution')
        end
    end
end


