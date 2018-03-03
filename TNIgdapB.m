function [ choi_ml_vec,solution, costs ] = gdapB( A,n )
%gdapB projected gradient descent with backtracking
%   Detailed explanation goes here
    d = sqrt(sqrt(size(A)));
    d = d(2);
    

%     n = reshape(n,[],2*d*d);% normalise clicks
%     for i=1:d*d
%         n(:,i) = n(:,i)/sum(n:,i)
%     end
    n = reshape(n,[],1);
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
    gamma = 0.01;%1e-6; % higher means more demanding line search. Boyd and Vandenberghe suggest between 0.01 and 0.3
    
%     Lscale = norm(gradient(A,n,choi_init));
    
    mu = 1; % inverse learning rate
    for i=1:1e10
%         mu = 1.05*mu;
%         i;
%         costs(i)     = 0; % just debugging
        costs(i)     = cost(A,n,solution{i}); % this not strictly necessary and quite expensive
%         costs(end)
        G = gradient(A,n,solution{i});
%         D{i}         = CPTP_project(solution{i}-(1/mu)*G, MdagM, Mdagb)-solution{i};
        D         = CPTNI_project(solution{i}-(1/mu)*G, MdagM, Mdagb)-solution{i};
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
        while cost(A,n,solution{i}+alpha*D) > cost(A,n,solution{i}) +  gamma*alpha*(D'*G)
            alpha = 0.8 * alpha ;
        end
%         solution{i+1} = solution{i} + alpha*D{i};
        solution{i+1} = solution{i} + alpha*D;
        if norm(solution{i+1}-solution{i})<1e-6
            break
        end
%         if cost(A,n,solution{i})-cost(A,n,solution{i+1})<1e-9
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
    plot(costs)
    hold on
    choi_ml_vec = CPTNI_project(solution{end}, MdagM, Mdagb);
end


