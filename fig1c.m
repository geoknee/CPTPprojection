d = 4;
A = PM_minimal(d);

% precompute matrices for TP_project
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

        
choi_ground     = randomCPTP_quasi_pure(d,0.9);
choi_ground_vec = reshape(choi_ground,[],1);

p               = real(A*choi_ground_vec);
%         p               = p/sum(p);
% n             = mnrnd(1e4,p)';
% n             = n/sum(n); % activate for multinomial noise

n               = p; % noiseless scenario
n               = n./sum(n); % noiseless scenario


if sum(n) == 1
    N = 10^12; % catch noiseless case this way
else       
    N = sum(n);
    n = n/N;
end
choi_init = sparse(eye(d*d)/d);
choi_init = reshape(choi_init,[],1);
solution  = {choi_init};
inside_cost(1) = cost(A,n,solution{1});
mu = 1000;%1.5/(d*d)
    for i=1:5

        G = gradient(A,n,solution{i});

        outside_point = solution{i}-(1/mu)*G;
        outside_cost(i) = cost(A,n,outside_point);
        inside_point  = CPTP_project(outside_point,MdagM, Mdagb);
        inside_cost(i+1) = cost(A,n,inside_point);
        solution{i+1} = inside_point;
        
%         D         = CPTP_project(solution{i}-(1/mu)*G, MdagM, Mdagb)-solution{i};
% 
%         alpha = 1;
%         new_cost = cost(A,n,solution{i});
%         B = new_cost + gamma*alpha*(D'*G);
%         while cost(A,n,solution{i}+alpha*D) > B  
%             alpha = 0.5 * alpha;  % less crude
%             B = new_cost + gamma*alpha*(D'*G);
%             if alpha < 1e-15
%                 break
%             end

        end

%         solution{i+1} = solution{i} + alpha*D;
% 
%         if (old_cost - new_cost)  < 1e-10
% 
%             break
%         end
%         old_cost = new_cost;

fig = figure;
set(fig,'position', [10 10 600 200])
ax1 = subplot(2,1,1);
scatter((1:6),inside_cost,'k');
yticks(ax1,[])
xticks(ax1,[])
xlim(ax1,[0,6])
ax1.XAxisLocation = 'origin';
ax1.YAxisLocation = 'origin';
ax2 = subplot(2,1,2);

scatter((1.5:5.5),outside_cost,'ks'); 
xticks(ax2,[1,2,3,4,5])
yticks(ax2,[])
xlim(ax2,[0,6])
saveas(gcf,'./plots/fig1c.eps','epsc')


