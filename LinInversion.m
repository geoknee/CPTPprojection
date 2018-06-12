function [ choi_LI_vec ] = LinInversion( A,n ) % non-interative so now convergence profile
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    d = sqrt(sqrt(size(A)));
    d = d(2);
    
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
    


%     A = full(A); % remove if using lsqminnorm
%     choi_LI_vec = pinv(A)*n;
%     choi_LI_vec = PSD_project(pinv(A)*n); % TODO check J distance ok for
%     not TP maps
%     choi_LI_vec = CPTP_project(pinv(A)*n, MdagM, Mdagb);
%     choi_LI_vec = lsqminnorm(A,n);
%     choi_LI_vec = PSD_project(A\n);
%     choi_LI_vec = A\n;
%     choi_LI     = reshape(choi_LI_vec,[],d*d);
%     choi_LI     = d*choi_LI/trace(choi_LI); % correct normalisation
%     choi_LI_vec     = reshape(choi_LI,[],1);
    choi_LI_vec = CPTP_project(A\n, MdagM, Mdagb);
    
%     choi_LI     = reshape(choi_LI_vec,[],d*d);
%     norm(partial_trace(choi_LI)-eye(d))
%     if prod(eig(choi_LI)>-1e-3)
%         fprintf('CP')
%     else
%         fprintf('not CP')
%         min(eig(choi_LI))
%     end
%         
%     if norm(partial_trace(choi_LI)-eye(d))<1e-4
%         fprintf('TP')
%     else
%         fprintf('not TP')
%         norm(partial_trace(choi_LI)-eye(d))        
%     end
    
end

