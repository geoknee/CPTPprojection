function [ choi_ml_vec, solution, costs  ] = DIA( A,n )
%UNTITLED Summary of this function goes here
%   here we follow Anis and Lvovsky NJP 2012
    d = sqrt(sqrt(size(A)));
    d = d(2);
    choi_init = eye(d*d)/d;
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
    e = 0.01;    
    for k=1:5e3
        costs(k)     = cost(A,n,solution{k});
        rho = reshape(solution{k},[],d*d) ;
        R   = reshape(gradient(A, n , solution{k}),[],d*d);
        R = (e*R)+(1-e)*eye(d*d); % dilute
        rho_new = R*rho*R;
        if prod(isfinite(partial_trace(rho_new))) == 0
            break
        end
        % TODO catch singular matrices, too
        lambda = sqrtm(partial_trace(rho_new));
        Lambda = kron(eye(d),lambda);
        if rcond(Lambda) < 1e-10
            break
        end
        LI = inv(Lambda);
        rho_new = LI*rho_new*LI;
%         rho_new = d*rho_new/trace(rho_new);
        solution{k+1} = reshape(rho_new,[],1);
        if k>20
            if var(costs(k-20:k)) < 1e-12
                break
            end
        end
    end
    choi_ml_vec = solution{end};
end

