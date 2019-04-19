function [ choi_ml_vec,solution, costs ] = gdapB( A,n )
%gdapB projected gradient descent with backtracking
    d = sqrt(sqrt(size(A)));
    d = d(2);
    
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
    
    % unbiased initialisation
    choi_init = sparse(eye(d*d)/d);        
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
    gamma = 0.3; % higher means more demanding line search. Boyd and Vandenberghe suggest between 0.01 and 0.3
    
    
    mu = 1.5/(d*d); % inverse learning rate
    old_cost = 1e10;
    costs = 0;
    cs{1} = 1e10;
    for i=1:1e10


        G = gradient(A,n,solution{i});
        D         = CPTP_project(solution{i}-(1/mu)*G, MdagM, Mdagb)-solution{i};
        alpha = 1;
        new_cost = cost(A,n,solution{i});
        B = new_cost + gamma*alpha*(D'*G);
        while cost(A,n,solution{i}+alpha*D) > B  
            alpha = 0.5 * alpha;  % less crude
            B = new_cost + gamma*alpha*(D'*G);
            if alpha < 1e-15
                break
            end
        end
        solution{i+1} = solution{i} + alpha*D;  
        if (old_cost - new_cost)  < 1e-10
            break
        end
        old_cost = new_cost;
    end
    choi_ml_vec = solution{end};
end


