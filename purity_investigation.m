clear;close all;
figure
for d=2:5
    d
    
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
    

    for m=1:d^2
        for i=1:1000

            purity_gaussian(i)=real(sum(eig(randomCPTP(d,m)).^2)/d^2);

        end
        mean_purity_gaussian(m) = mean(purity_gaussian);
        std_purity_gaussian(m)  = std(purity_gaussian);
    end
    for i=1:1000

        choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
        choi_ground_vec = reshape(choi_ground,[],1);
        choi_ground_vec = CPTP_project(choi_ground_vec, MdagM, Mdagb);
        choi_ground     = reshape(choi_ground_vec,[],d*d);

        rank_square(i)  = rank(choi_ground);
        purity_square(i)= sum(eig(choi_ground).^2/d^2);
        
        
        choi_quasi      = randomCPTP_quasi_pure(d,0.9);
        
        rank_quasi(i)  = rank(choi_quasi);
        purity_quasi(i)= sum(eig(choi_quasi).^2/d^2);
        
    end
    mean_rank_square = mean(rank_square);
    std_rank_square  = std(rank_square);

    mean_purity_square = mean(purity_square);
    std_purity_square = std(purity_square);
    
    mean_rank_quasi = mean(rank_quasi);
    std_rank_quasi  = std(rank_quasi);
    
    mean_purity_quasi = mean(purity_quasi);
    std_purity_quasi = std(purity_quasi);
    
    ax = gca;
    ax.ColorOrderIndex = d-1;

    errorbar(1:d^2,mean_purity_gaussian,std_purity_gaussian,'LineWidth',2,'DisplayName',num2str(d)); hold on
    
    ax = gca;
    ax.ColorOrderIndex = d-1;
    errorbar(mean_rank_square,mean_purity_square,std_purity_square,std_purity_square,std_rank_square,std_rank_square,'d','LineWidth',2,'DisplayName',num2str(d))
    errorbar(mean_rank_quasi,mean_purity_quasi,std_purity_quasi,std_purity_quasi,std_rank_quasi,std_rank_quasi,'d','LineWidth',2,'DisplayName',num2str(d))

    end
%     print(gcf,'-depsc2','-loose',['./plots/puritygaussiand',num2str(d),'.eps'])
xlabel('rank')
ylabel('purity')
legend('show','Location','northeast')
set(gca,'fontsize',18)
ylim([0,1])
% set(gca,'yscale','log')
saveas(gcf,['./plots/puritygaussiand',num2str(d),'.eps'],'epsc')
    
    
% figure
%     
%     
% for d=2:5
%     d
%        
%     M = zeros([d*d,d*d*d*d]);
%     for i=1:d
%         e = zeros(1,d);
%         e(i)  = 1;
%         B = kron(eye(d),e); 
%         M = M + kron(B,B);
%     end
%     MdagM = sparse(M'*M);
%     b = sparse(reshape(eye(d),[],1));
%     Mdagb = sparse(M'*b);
%         
%     for i=1:100
% 
%         choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
%         choi_ground_vec = reshape(choi_ground,[],1);
%         choi_ground_vec = CPTP_project(choi_ground_vec, MdagM, Mdagb);
%         choi_ground     = reshape(choi_ground_vec,[],d*d);
%         
% %         purity_square(i)=real(sum(eig(choi_ground).^2)/d^2)
%     end
%     mean_rank_square(d) = mean(rank(choi_ground));
%     std_rank_square(d)  = std(rank(choi_ground));
%     
%     mean_purity_square(d) = mean(sum(eig(choi_ground).^2/d^2));
%     std_purity_square(d) = std(sum(eig(choi_ground).^2/d^2));
% 
% %     plot(purity_square); hold on
% end
% 
% errorbar(mean_rank_square,mean_purity_square,std_rank_square,std_rank_square,std_purity_square,std_purity_square,'d')
% xlabel('rank')
% ylabel('purity')
% legend('d=2','d=3','d=4','d=5')
% set(gca,'fontsize',18)
% saveas(gcf,['./plots/puritysquare',num2str(d),'.eps'],'epsc')
% % print(gcf,'-depsc2','-loose',['./plots/puritysquare',num2str(d),'.eps'])
% 
% 
% 
