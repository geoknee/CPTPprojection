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
    


    A = full(A); % remove if using lsqminnorm
%     choi_LI_vec = pinv(A)*n;
%     choi_LI_vec = PSD_project(pinv(A)*n); % TODO check J distance ok for
%     not TP maps
%     choi_LI_vec = CPTP_project(pinv(A)*n, MdagM, Mdagb);
%     choi_LI_vec = lsqminnorm(A,n);
    choi_LI_vec = PSD_project(A\n);
%     choi_LI_vec = CPTP_project(A\n, MdagM, Mdagb);

    
end

