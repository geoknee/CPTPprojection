function [ choi_ml_vec, solution, costs  ] = DIA( A,n )
%DIA: diluted RrhoR algorithm 
%   here we follow Anis and Lvovsky NJP 2012
% also Fedorov et al NJP 2014

    d = sqrt(sqrt(size(A)));
    d = d(2);
    choi_init = eye(d*d)/d;
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
    costs = zeros(1e4);
    for k=1:1e4

%         k
%         costs(k)
        rho = reshape(solution{k},[],d*d) ;
        rho = 0.5*rho + 0.5*rho'; %enforce hermitian
%         R = 0;
%         for i=1:length(n)
%             Rn   = conj(reshape(A(i,:),[],d*d)); % conj is important
%             Rd   = real(trace(Rn*rho));
%             R = R + n(i)*Rn*1.0/Rd; %n(i) missing from Anis paper since they do not bin their data, all outcomes unique.
%             % weirdly the minus sign makes a lot of difference...
%         end
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
%             Y  = kron(partial_trace(0.5*R*rho+0.5*rho*R),eye(d));
%             LI = (1+e)*eye(d*d)- e * Y; % equation (A6) of Anis et al
        end

        rho_new = LI*Rp*rho*Rp*LI;
        solution{k+1} = reshape(rho_new,[],1);
%         norm(rho_new-rho)

        if norm(rho_new-rho,'fro')<1e-6 %warning, making this very strict will result in very large datafiles
%             k
%             costs(end)
%             partial_trace(rho_new)
%             eig(rho_new)
            break
        end
        costs(k+1)     = cost(A,n,solution{k+1});
%         if abs(costs(k)-costs(k+1))<1e-9
%             break
%         end
        
    end
    choi_ml_vec = solution{end};
end

