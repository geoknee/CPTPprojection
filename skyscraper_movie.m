

d = 4;

% choi_ground     = randomCPTP_quasi_pure(d,0.9);
choi_ground     = randomCPTP(d,d*d);
choi_ground_vec = reshape(choi_ground,[],1);

A = PM_minimal(d);

p               = real(A*choi_ground_vec);
n               = p; % noiseless scenario

[choi_ml_vec, solution, costs] = gdapB(A,n);

v = VideoWriter('./plots/movies/FULLRANKss.avi');
open(v)
for o=1:floor(length(solution)/3)
    choi_ml_vec = solution{o}
    choi_ml = reshape(choi_ml_vec,[],d*d);

    error = trace_dist(choi_ground/trace(choi_ground),choi_ml/trace(choi_ml))

    figure('Position',[1 0 1300 1300]);
    subplot(2,2,1)
    
    b1 = bar3(real(choi_ground));
    title('Ground truth (Real)')

    for k = 1:length(b1)
        zdata = b1(k).ZData;
        b1(k).CData = zdata;
        b1(k).FaceColor = 'interp';
    end

    xlabel('i')
    ylabel('j')
    xticks([1,d*d/2,d*d])
    yticks([1,d*d/2,d*d])

    ylim([0,d*d+0.5])
    xlim([0,d*d+0.5])
    zlim([-1,1])
    zlabel('Re$C_{\mathcal{E}}$','Interpreter','latex')
    set(gca,'fontsize',18)
    
    
    subplot(2,2,3)
    

    b3 = bar3(real(choi_ml));
    title(['Estimate (Real) k = ',num2str(o)])

    for k = 1:length(b3)
        zdata = b3(k).ZData;
        b3(k).CData = zdata;
        b3(k).FaceColor = 'interp';
    end

    xlabel('i')
    ylabel('j')
    xticks([1,d*d/2,d*d])
    yticks([1,d*d/2,d*d])

    ylim([0,d*d+0.5])
    xlim([0,d*d+0.5])
    zlim([-1,1])
    zlabel('Re$C_{\mathcal{E}}$','Interpreter','latex')
    set(gca,'fontsize',18)
    
        
    subplot(2,2,2)
    b2 = bar3(imag(choi_ground));
    title('Ground truth (Imag)')

    for k = 1:length(b2)
        zdata = b2(k).ZData;
        b2(k).CData = zdata;
        b2(k).FaceColor = 'interp';
    end

    xlabel('i')
    ylabel('j')
    xticks([1,d*d/2,d*d])
    yticks([1,d*d/2,d*d])
    ylim([0,d*d+0.5])
    xlim([0,d*d+0.5])
    zlim([-1,1])
    zlabel('Im$C_{\mathcal{E}}$','Interpreter','latex')

    set(gca,'fontsize',18)
    
    subplot(2,2,4)
    b4 = bar3(imag(choi_ml));
    title(['Estimate (Imag) error = ',num2str(error)])

    for k = 1:length(b4)
        zdata = b4(k).ZData;
        b4(k).CData = zdata;
        b4(k).FaceColor = 'interp';
    end

    xlabel('i')
    ylabel('j')
    xticks([1,d*d/2,d*d])
    yticks([1,d*d/2,d*d])
    ylim([0,d*d+0.5])
    xlim([0,d*d+0.5])
    zlim([-1,1])
    zlabel('Im$C_{\mathcal{E}}$','Interpreter','latex')

    set(gca,'fontsize',18)
    
    F=getframe(gcf);
    writeVideo(v,F)
   
    writeVideo(v,F) % write same frame 2 times
    
%     saveas(gcf,'./plots/movies/ss.png')
    close
    % saveas(gcf,'./plots/ss.eps','epsc')
end

close(v)

