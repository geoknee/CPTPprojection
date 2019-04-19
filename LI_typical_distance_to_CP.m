dmax = 4;
Npowmin = 3;
Npowmax = 6;
ensemble_size = 1000;
for d = 4:4
    A = PM_minimal(d);
    

    for l=1:ensemble_size
            fprintf(char(10));
            fprintf('%d ', l); 

            choi_ground     = randomCPTP_quasi_pure(d,0.9);
            choi_ground_vec = reshape(choi_ground,[],1);

            p               = real(A*choi_ground_vec);

    %         for Npow=[1,2,3,4,5,6,7,8,Inf] % above Npow=9 the memory requirements are huge for simulating multinomial noise
            for Npow=Npowmin:Npowmax                            
                N = 10^Npow;

                if isinf(N)
                    n           = p;
                else
                    pmat           = reshape(p,[],2*d*d); % need an object with n_measurement_outcomes columns
                    pmat           = pmat./sum(pmat,2); % does not look necessary but it is useful to avoid near misses where probs sum to 1-e.
                    nmat           = mnrnd(N,pmat);
                    nmat        = nmat./sum(nmat,2); % proper normalisation so that sum(nmat,2)=1
                    n           = reshape(nmat,[],1);
                end
                
                choi_li = LinInversion(A,n);
                choi_li_matrix = reshape(choi_li,[],d*d);
                choi_li_matrix = (choi_li_matrix'+choi_li_matrix)/2;
                smallest_ev(d,Npow,l)           = min(eig(choi_li_matrix));
                frob_dist_to_TP(d,Npow,l)       = norm(partial_trace(choi_li_matrix)-eye(d),'fro');       
            end
    end
end
figure('Position',[1 0 800 250]);
hold on;
m  = mean(smallest_ev,3);
sd = std(smallest_ev,0,3);
for d=2:dmax
    errorbar(Npowmin:Npowmax,m(d,Npowmin:Npowmax),sd(d,Npowmin:Npowmax),'LineWidth',2)
end
xlim([Npowmin-0.5,Npowmax+0.5])
ylabel('smallest eigenvalue')
xlabel('log(N)')
xticks(Npowmin:Npowmax)
set(gca,'YScale','log');
figure('Position',[1 0 800 250]);

hold on;
for Npow=Npowmin:Npowmax
    subplot(1,2,1)
    hold on;box on;
    xlabel('smallest eigenvalue')
    ylabel('number')
    h = histogram(real(smallest_ev(d,Npow,:)),'DisplayName',num2str(10^Npow));
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    h.EdgeAlpha = 1;
    subplot(1,2,2)
    xlabel('distance to TP')
    ylabel('number')
    hold on;
    hold on; box on;
    h = histogram(real(frob_dist_to_TP(d,Npow,:)),'DisplayName',num2str(10^Npow));
    h.EdgeColor = 'none';
    h.FaceAlpha = 1;
    h.EdgeAlpha = 1;
end
lgd = legend('show');
title(lgd,'N')
yl = ylim;
plot([0,0],yl,'-k')

saveas(gcf,['./plots/d',num2str(d),'smallestEV.png'])
saveas(gcf,['./plots/d',num2str(d),'smallestEV.eps'],'epsc')
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2','-loose',['./plots/ALTd',num2str(d),'smallestEV.eps'])
    