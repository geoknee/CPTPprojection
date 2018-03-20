function [ g ] = gradient( A, n ,choi_vec)
%gradient: calulates the gradient of the cost function we try to minimize
% A             : (n_measurements x d^4 ) dimensional feature matrix
% n             : observed counts n
% choi_vec      : (d^4 x 1) vectorised choi matrix (candidate process)
% g             : (d^4 x 1) is the gradient of the negative log-likelihood 
  
    p   = real(A*choi_vec);
    
%     n   = n./sum(n);
%     p   = p./sum(p);
    
    g = -A'*(n./p);

end
