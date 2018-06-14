clear all;
dmax = 4;
Npowmin = 2;
Npowmax = 5;
ensemble_size = 10;
drange = 2:7;
parpool(4)
for d = drange
    fprintf(char(10));
    fprintf('%d ', d); 
    A = PM_minimal(d);
    for l=1:ensemble_size
            fprintf(char(10));
            fprintf('%d ', l); 
            % generate random ground truth
            choi_ground     = randomCPTP(d,d*d);
%             choi_ground     = randomCPTP_quasi_pure(d,0.9);
            choi_ground_vec = reshape(choi_ground,[],1);
            p               = real(A*choi_ground_vec);

            parfor Npow=Npowmin:Npowmax                            
                N = 10^Npow;

                if isinf(N)
                    n           = p;
                else
                    pmat           = reshape(p,[],2*d*d); % need an object with n_measurement_outcomes columns
                    pmat           = pmat./sum(pmat,2); % does not look necessary but it is useful to avoid near misses where probs sum to 1-e.
                    nmat           = mnrnd(N,pmat);
                    nmat           = nmat./sum(nmat,2); % proper normalisation so that sum(nmat,2)=1
                    n              = reshape(nmat,[],1);
                end
                cvx_solver mosek
                choi_cvx = cvx_wrapper(A,n);
                if isnan(choi_cvx)
                    fprintf('fail    ')
                    fail(d,Npow,l)=1;
                else
                    fprintf('success ')
                    fail(d,Npow,l)=0;
                end
       
            end
    end
end

%%
fail_freq = sum(fail,3)./size(fail,3);
% bar3(fail_freq(drange,Npowmin:Npowmax))
bar3(drange,fail_freq(drange,:))
ylabel('d')
xlabel('logN')
zlabel('failure frequency')

saveas(gcf,['./plots/cvx_failure.png'])
saveas(gcf,['./plots/cvx_failure.eps'],'epsc')
savefig(gcf,['./plots/cvx_failure'])