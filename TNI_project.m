function [ projected_choi_vec ] = TNI_project( choi_vec)
%TP_projectt: project a matrix into the set of trace nonincreasing maps
% choi_vec          : is a vector with dimensions (d^4 x 1).
%projected_choi_vec : (d^4 x 1) which represents a vectorised TNI Choi
%                   : matrix
    d = sqrt(sqrt(size(choi_vec)));
    d = d(1);
    
    choi = reshape(choi_vec,[],d*d);
    
    Y = partial_trace(choi);
    
    
   
    [V,D] = eig(eye(d)-Y);
    D = max(real(D),0);
    Y_proj = eye(d) - V*D*V';
    
    
    choi_proj = choi + kron(eye(d),Y_proj-Y);
   
    projected_choi_vec = reshape(choi_proj,[],1);
    
end
