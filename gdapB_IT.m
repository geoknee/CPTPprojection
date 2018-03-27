function [ choi_ml_vec,solution, costs ] = gdapB( A,n )
%gdapB projected gradient descent with backtracking
%   Detailed explanation goes here
    d = sqrt(sqrt(size(A)));
    d = d(2);
    

%     n = reshape(n,[],2*d*d);% normalise clicks
%     for i=1:d*d
%         n(:,i) = n(:,i)/sum(n:,i)
%     end
%     n = reshape(n,[],1);
    % precompute matrices for TP_project
    M = zeros([d*d,d*d*d*d]);
    for i=1:d
        e = zeros(1,d);
        e(i)  = 1;
        B = kron(speye(d),e); 
        M = M + kron(B,B);
    end
    MdagM = sparse(M'*M);
    b = reshape(speye(d),[],1);
    Mdagb = sparse(M'*b);
    
    choi_init = sparse(eye(d*d)/d);
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
%     stepsize      = 1.0/(1e3*d);
    gamma = 0.3;%1e-6; % higher means more demanding line search. Boyd and Vandenberghe suggest between 0.01 and 0.3
    
%     Lscale = norm(gradient(A,n,choi_init));
    
    mu = 1.5/(d*d); % inverse learning rate
    old_cost = 1e10;
    costs = 0;
    q = 1e-3;
    n = (1-q)*n + q/d;
    for i=1:1e10
%         mu = 1.05*mu;
%         i;
        solution{i} = (1-q)*solution{i}+q*reshape(eye(d*d),[],1)/d;
        G = gradient(A,n,solution{i});
%         norm(G)
%         if norm(G)<1e-6
%             break
%         end
        
%         % glancy stopping criterion NJP 14 095017
%         
%         r(i) = max(eig(reshape(-G,[],d*d))) % TODO glancy has -N term. we are omitting
%         if r(i)<1e-6*0.5*(d*d-1)
%             break
%         end
%         
        
%         D{i}         = CPTP_project(solution{i}-(1/mu)*G, MdagM, Mdagb)-solution{i};
        D         = CPTP_project(solution{i}-(1/mu)*G, MdagM, Mdagb)-solution{i};
%         sum(svd(D{i}))
%         if sum(svd(D{i}))<1e-15 % no point using trace norm because these
%         are vectors
%         norm(D{i})
%         if norm(D{i})<(1e-4)%*1/d
% %          if sum(svd(D{i}))<1e-4
% %             i
% %             costs(end)
%             break
%         end
%         if norm(D{i})<1e-10   
%             break
%         end
        alpha = 1;
        new_cost = cost(A,n,solution{i});
        B = new_cost + gamma*alpha*(D'*G);
        while cost(A,n,solution{i}+alpha*D) > B  
            alpha = 0.5 * alpha;  % less crude
            B = new_cost + gamma*alpha*(D'*G);
            if alpha < 1e-15
                break
            end
%             alpha = 0.2 * alpha ; % more crude
        end
%         solution{i+1} = solution{i} + alpha*D{i};
        solution{i+1} = solution{i} + alpha*D;
        solution{i+1} = (solution{i+1} - (q/d)*reshape(eye(d*d),[],1))/(1-q);
% 
%         if norm(solution{i+1}-solution{i})<5e-5 % criterion in solution space rather than costs seems to work better for pgd
%             break
%         end
%         
        if (old_cost - new_cost)  < 1e-10
%             old_cost-new_cost
%         if (new_cost)/old_cost > 1- 1e-5
%             new_cost
            break
        end
        old_cost = new_cost;
%         costs(i)     = 0; % just debugging
%         costs(i+1)     = cost(A,n,solution{i+1}); % this not strictly necessary and quite expensive
%         costs(end)
%         if abs(costs(i)-costs(i+1))<1e-9
%             break
%         end
%         if i>1
% %             if var(costs(i-10:i)) < 1e-13
%             if norm(alpha*D{i})<1e-14
%                 break
%             end
%         end
%     costs(end)
%         if i>1
%             if norm(solution{i}-solution{i-1})<1e-12
% %                 i
%                 break
%             end
%         end
        
        
    end
%     plot(costs)
%     hold on
    choi_ml_vec = CPTP_project(solution{end}, MdagM, Mdagb);
%     choi_ml_vec = solution{end};
end

