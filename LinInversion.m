function [ choi_LI_vec ] = LinInversion( A,n ) % non-interative so no convergence profile
%Linear Inversion with a Final Projection (called LIFP in the paper)
%   Uses matlab backslash to solve unconstrained, and then uses a single CPTP projection
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
    
    choi_LI_vec = CPTP_project(A\n, MdagM, Mdagb);
    
end

