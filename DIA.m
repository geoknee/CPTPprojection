function [ choi_ml_vec, solution, costs  ] = DIA( A,n )
%DIA: diluted RrhoR algorithm 
%   here we follow Anis and Lvovsky NJP 2012
% also Fedorov et al NJP 2014

    d = sqrt(sqrt(size(A)));
    d = d(2);
    choi_init = eye(d*d)/d;
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};

    for k=1:2000
        costs(k)     = cost(A,n,solution{k});
%         k
%         costs(k)
        rho = reshape(solution{k},[],d*d) ;
        R = 0;
        for i=1:length(n)
            Rn   = reshape(A(i,:),[],d*d);
            Rd   = real(trace(Rn*rho));
            R = R + n(i)*Rn*1.0/Rd; %n(i) missing from Anis paper since they do not bin their data, all outcomes unique
        end
        e = 1; 
        R = (e*R)+(1-e)*eye(d*d); % dilute
        lambda  =  sqrtm(partial_trace(R*rho*R));
        LI      =  kron(inv(lambda),eye(d));
        while cost(A,n,reshape(LI*R*rho*R*LI,[],1)) > cost(A,n,reshape(rho,[],1))
            e       =  e * 0.2;
            if e < 1e-100
                break
            end
            R       =  (e*R)+(1-e)*eye(d*d);
            lambda  =  sqrtm(partial_trace(R*rho*R));
            lambda  = eye(d);
            LI      =  kron(inv(lambda),eye(d)); 
%             Y  = kron(partial_trace(0.5*R*rho+0.5*rho*R),eye(d));
%             LI = (1+e)*eye(d*d)- e * Y; % equation (A6) of Anis et al
        end

        rho_new = LI*R*rho*R*LI;

        solution{k+1} = reshape(rho_new,[],1);
%         norm(rho_new-rho)

        if norm(rho_new-rho)<1e-9
%             k
%             costs(end)
%             partial_trace(rho_new)
%             eig(rho_new)
            break
            
        end
        
    end
    choi_ml_vec = solution{end};
end

