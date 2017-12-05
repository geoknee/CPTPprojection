
function [ choi_ml_vec,solution, costs ] = gdapB( A,n )
% with backtracking
%   Detailed explanation goes here
    d = sqrt(sqrt(size(A)));
    d = d(2);
    choi_init = eye(d*d)/d;
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
%     stepsize      = 1.0/(1e3*d);
    gamma = 1e-7;
    for i=1:5e3
        i;
        costs(i)     = cost(A,n,solution{i});
        D{i}  = CPTP_project(solution{i}-gradient(A,n,solution{i}))-solution{i};
%         if norm(D{i})<1e-10
%             break
%         end
        alpha = 1;
        while cost(A,n,solution{i}+alpha*D{i}) > cost(A,n,solution{i}) +  gamma*alpha*(D{i}'*gradient(A,n,solution{i}))
            alpha = 0.2 * alpha ;
        end
        solution{i+1} = solution{i} + alpha*D{i};
        if i>20
            if var(costs(i-20:i)) < 1e-5
                break
            end
        end
    end

    choi_ml_vec = solution{end};
end

