function [ projected_choi_vec ] = TNI_project( choi_vec)
%TP_projectt: project a matrix into the set of trace nonincreasing maps
% choi_vec          : is a vector with dimensions (d^4 x 1).
%projected_choi_vec : (d^4 x 1) which represents a vectorised TNI Choi
%                   : matrix
    d = sqrt(sqrt(size(choi_vec)));
    d = d(1);
    
    choi = reshape(choi_vec,[],d*d);
    
    choi = (choi+choi')/2;
    
    Y = partial_trace(choi);
    
%     Y = (Y+Y')/2;
%     [V,D] = eig(Y,'nobalance');
    [V,D] = eig(Y);
    D = min(real(D),1);
    Y_proj = V*D*V';
    
    
    choi_proj = choi + kron(Y_proj-Y,eye(d))/d;
   
    projected_choi_vec = reshape(choi_proj,[],1);
    
end
