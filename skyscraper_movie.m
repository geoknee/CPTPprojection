

d = 4;

choi_ground     = randomCPTP_quasi_pure(d,0.9);
choi_ground_vec = reshape(choi_ground,[],1);

A = PM_minimal(d);

p               = real(A*choi_ground_vec);
n               = p; % noiseless scenario

[choi_ml_vec, solution, costs] = gdapB(A,n);

v = VideoWriter('./plots/movies/ss.avi');
open(v)
for o=1:floor(length(solution)/5)
    choi_ml_vec = solution{o}
    choi_ml = reshape(choi_ml_vec,[],d*d);

    trace_dist(choi_ground/trace(choi_ground),choi_ml/trace(choi_ml))

    figure('Position',[1 0 1300 500]);
    subplot(1,2,1)

    b1 = bar3(real(choi_ml));

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
    subplot(1,2,2)
    b2 = bar3(imag(choi_ml));

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
    
    F=getframe(gcf);
    writeVideo(v,F)
   
    writeVideo(v,F) % write same frame 2 times
    
%     saveas(gcf,'./plots/movies/ss.png')
    close
    % saveas(gcf,'./plots/ss.eps','epsc')
end

close(v)

