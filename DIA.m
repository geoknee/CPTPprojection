function [ choi_ml_vec, solution, costs  ] = DIA( A,n )
%DIA: diluted RrhoR algorithm 
%   here we follow Anis and Lvovsky NJP 2012
% also Fedorov et al NJP 2014

    d = sqrt(sqrt(size(A)));
    d = d(2);
    
    % unbiased initialisation
    choi_init = sparse(eye(d*d)/d);
    
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
    costs=0;
    old_cost = 1e10;
    for k=1:1e4

        rho = reshape(solution{k},[],d*d) ;
        rho = 0.5*rho + 0.5*rho'; %enforce hermitian
        R = reshape(gradient(A,n,solution{k}),[],d*d);
        e = 1; 
        Rp = (e*R)+(1-e)*eye(d*d); % dilute
        lambda  =  sqrtm(partial_trace(Rp*rho*Rp));
        LI      =  kron(inv(lambda),eye(d));
        while cost(A,n,reshape(LI*Rp*rho*Rp*LI,[],1)) > cost(A,n,reshape(rho,[],1))
            e       =  e * 0.5;
            if e < 1e-100
                break
            end
            Rp       =  (e*R)+(1-e)*eye(d*d);
            lambda  =  sqrtm(partial_trace(Rp*rho*Rp));
            LI      =  kron(inv(lambda),eye(d)); 
        end

        rho_new = LI*Rp*rho*Rp*LI;
        solution{k+1} = reshape(rho_new,[],1);
        new_cost = cost(A,n,solution{k+1});
        
        p = real(A*solution{k+1});
        
        if (old_cost - new_cost)  < 1e-10
            break
        end
        old_cost = new_cost;
        
    end
    choi_ml_vec = solution{end};
end

