function [ choi_ml_vec,solution, costs ] = gdapM( A,n )
%gdapM projected gradient descent with momentum
    d = sqrt(sqrt(size(A)));
    d = d(2);
    

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



    gamma = 1/(2*d);
    xi    = 0.95; % inertia
    PM    = ceil(log(cost(A,n,choi_init)));
    T{1}  = 0; % momentum
    
    for k=1:1e10
        solution{k}  = CPTP_project(solution{k},MdagM,Mdagb);
        costs(k)     = cost(A,n,solution{k}); % this not strictly necessary and quite expensive
        G = gradient(A,n,solution{k});
        CM = ceil(log(cost(A,n,solution{k})));
        if CM < PM
            xi = 1-(1-xi)*0.95; % update inertia
            PM = CM;
        end

        T{k+1}        = xi*T{k}-gamma*G;
        solution{k+1} = solution{k}+T{k+1};
%         
        if k>20
            if var(costs(end-20:end))<1e-12
                break
            end
        end

    end
    plot(costs)
    hold on;
    choi_ml_vec = CPTP_project(solution{end},MdagM,Mdagb);


