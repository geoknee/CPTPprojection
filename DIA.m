function [ choi_ml_vec, solution, costs  ] = DIA( A,n )
%DIA: diluted RrhoR algorithm 
%   here we follow Anis and Lvovsky NJP 2012
% also Fedorov et al NJP 2014

    d = sqrt(sqrt(size(A)));
    d = d(2);
    choi_init = eye(d*d)/d;
    choi_init = reshape(choi_init,[],1);
    solution  = {choi_init};
    
%         % precompute matrices for TP_project
%     M = zeros([d*d,d*d*d*d]);
%     for i=1:d
%         e = zeros(1,d);
%         e(i)  = 1;
%         B = kron(eye(d),e); 
%         M = M + kron(B,B);
%     end
%     MdagM = sparse(M'*M);
%     b = sparse(reshape(eye(d),[],1));
%     Mdagb = sparse(M'*b);
%     
    
    e = 5e-4; % perhaps make this dimension dependent?
%     costs = zeros(1,1e2);
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
%         R   = reshape(gradient(A, n , solution{k}),[],d*d);
        R = (e*R)+(1-e)*eye(d*d); % dilute
%         R = (e*R)+eye(d*d); % dilute
        rho_new = R*rho*R;
%         if prod(isfinite(partial_trace(rho_new))) == 0
%             break
%         end
        % TODO catch singular matrices, too
        lambda = sqrtm(partial_trace(rho_new));
%         lambda = eye(d); % debugging, turn this part of algo off
%         Lambda = kron(eye(d),lambda);
        Lambda = kron(lambda,eye(d)); % DEBUG THIS WORKS
        % try hybridizing with my TP_project subroutine instead?
%         if rcond(Lambda) < 1e-10
%             break
%         end
        LI = inv(Lambda);
        rho_new = LI*rho_new*LI;
%         rho_new = d*rho_new / trace(rho_new);
%         partial_trace(rho_new)
%         rho_new = d*rho_new/trace(rho_new);
        solution{k+1} = reshape(rho_new,[],1);
%         solution{k+1} = TP_project(solution{k+1}, MdagM, Mdagb); % TRYING THIS OUT
%         costs(end)
%         norm(rho_new-rho)
        if norm(rho_new-rho)<1e-4
            break
        end
        
%         if k>20
%             if var(costs(k-20:k)) < 1e-12
%                 break
%             end
%         end
    end
    choi_ml_vec = solution{end};
end

