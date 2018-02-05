
function [ choi_ml_vec,solution, costs ] = gdapB( A,n )
% with backtracking
%   Detailed explanation goes here
    d = sqrt(sqrt(size(A)));
    d = d(2);
    
    % precompute matrices for TP_project
    M = zeros([d*d,d*d*d*d]);
    for i=1:d
        e = zeros(1,d);
        e(i)  = 1;
        B = kron(eye(d),e); 
        M = M + kron(B,B);
    end
    M = sparse(M);
    MdagM = M'*M; 
    
    
    choi_init = sparse(eye(d*d)/d);
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
%     stepsize      = 1.0/(1e3*d);
    gamma = 1e-7;
    mu = 1; % inverse learning rate?
    for i=1:5e3
%         i;
        costs(i)     = cost(A,n,solution{i});
        D{i}         = CPTP_project(solution{i}-(1/mu)*gradient(A,n,solution{i}), MdagM, M)-solution{i};
%         sum(svd(D{i}))
        if sum(svd(D{i}))<1e-5
%         if norm(D{i})<1e-3
% %             i
            break
        end
%         if norm(D{i})<1e-10   
%             break
%         end
        alpha = 1;
        while cost(A,n,solution{i}+alpha*D{i}) > cost(A,n,solution{i}) +  gamma*alpha*(D{i}'*gradient(A,n,solution{i}))
            alpha = 0.5 * alpha ;
        end
        solution{i+1} = solution{i} + alpha*D{i};
%         if i>1
% %             if var(costs(i-10:i)) < 1e-13
%             if norm(alpha*D{i})<1e-14
%                 break
%             end
%         end
    end

    choi_ml_vec = CPTP_project(solution{end}, MdagM, M);
end


