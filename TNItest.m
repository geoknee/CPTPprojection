% test my TNI projection. 
d = 2;
cvx_wins = 0;
kmax = 100;
for k=1:kmax
% k
choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
choi_ground_vec = reshape(choi_ground,[],1);
choi_ground_vec = PSD_project(choi_ground_vec);

choi_ground     = reshape(choi_ground_vec,[],d*d);

% eig(partial_trace(choi_ground))

cvx_solver mosek
% cvx_precision best
cvx_begin quiet

variable cvx_solution(d*d,d*d) complex

minimize(norm(cvx_solution-choi_ground,'fro'))

subject to 

eye(d) - partial_trace(cvx_solution) == semidefinite(d,d);

cvx_end

cvx_solution;

GK_solution     = reshape(TNI_project(choi_ground_vec),[],d*d);

discrepancies(k) = norm(GK_solution-cvx_solution,'fro');
cvx_distance = norm(choi_ground-cvx_solution,'fro');
 gk_distance = norm(GK_solution-choi_ground,'fro');
 if cvx_distance < gk_distance
     cvx_wins = cvx_wins + 1;
 end
%  cvx_win_percentage = cvx_wins/k
end

figure;
h1 = histogram(discrepancies); hold on;
% h2 = histogram(cvx_distances);
% h3 = histogram(gk_distances);
% legend('discrepancy','cvx dist','TNI dist')
xlabel('discrepancy between cvx and TNI project');
ylabel('counts')
mean(discrepancies)
