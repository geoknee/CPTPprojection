function [ c ] = cost( A, n , choi_vec )
%cost: calulates the cost function we try to minimize
% A             : (n_measurements x d^4 ) dimensional feature matrix
% n             : observed counts n
% choi_vec      : (d^4 x 1) vectorised choi matrix (candidate process)
% c             : (scalar) is the negative log-likelihood 

    p  = real(A*choi_vec); % the forward model
%     n  = n./sum(n);
%     p  = p./sum(p);


%     eps = 1e-16;
%     
%     if nnz(p)<length(p)
%         sprintf('p=0 occured')
%     end
%     
%     p(find(p<eps)) = eps;
    
    
    c = real(-(n.'*reallog(p)));
end

