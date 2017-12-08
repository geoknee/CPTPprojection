
function [ choi_ml_vec,solution, costs ] = gdapB( A,n )
% with backtracking
%   Detailed explanation goes here
    d = sqrt(sqrt(size(A)));
    d = d(2);
    choi_init = eye(d*d)/d;
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
%     stepsize      = 1.0/(1e3*d);
    gamma = 1e-7; % larger means slower
    costs = zeros(1,1e2);
%     D = cell(1,1e3);
    for i=1:1e2
        costs(i)  = cost(A,n,solution{i});
        D         = CPTP_project(solution{i}-gradient(A,n,solution{i}))-solution{i};
%         if norm(D{i})<1e-10
%             break
%         end
        alpha = 1;
        while cost(A,n,solution{i}+alpha*D) >= cost(A,n,solution{i}) +  gamma*alpha*(D'*gradient(A,n,solution{i}))
            alpha = 0.2 * alpha ;
        end
        solution{i+1} = solution{i} + alpha*D;
        if i>20
            if var(costs(i-20:i)) < 1e-5
                break
            end
        end
    end

%     choi_ml_vec = CPTP_project(solution{end});
    choi_ml_vec = solution{end};
end

