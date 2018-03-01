function [ choi_ml_vec] = cvx_wrapper( A,n )
    d = sqrt(sqrt(size(A)));
    d = d(2);
    
    M = zeros([d*d,d*d*d*d]);
    for i=1:d
        e = zeros(1,d);
        e(i)  = 1;
        B = kron(eye(d),e);
        M = M + kron(B,B);
    end
    b = reshape(eye(d),[],1);
    
    
%     cvx_solver sedumi
    cvx_begin %quiet
    variable cvx_choi(d*d,d*d) hermitian semidefinite
    variable P(d,d) hermitian semidefinite
    % todo TP constraint
    choi_vec_cvx = reshape(cvx_choi,d*d*d*d,[]);
    n = real(n);
    p_cvx        = real(A*choi_vec_cvx);
    maximize(n'*log(p_cvx))
    subject to
    % eye(d)-reshape(M*choi_vec_cvx,d,d) ==  P; % this is for TNI
    M*choi_vec_cvx == b; % this is for TP
    cvx_end
    cvx_time = toc;

    choi_ml_vec = choi_vec_cvx;
end

