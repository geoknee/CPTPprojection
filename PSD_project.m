function [ projected_choi_vec ] = PSD_project( choi_vec )
%PSD_projectt: project a matrix into the set of completely positive maps
% we rely on Choi's theorem, that the Choi matrix's positivity <--> the
% complete positivity of the map itself. Hence we use a projection of the
% Choi matrix onto the set of positive semidefinite matrices.
% choi_vec          : is a vector with dimensions (d^4 x 1).
%projected_choi_vec : (d^4 x 1) which represents a vectorised, positive Choi
%                   : matrix

    d = sqrt(size(choi_vec));
    d = d(1);
    choi = reshape(choi_vec,[],d);
       
    choi = (choi + choi')/2; % project onto Hermitian matrices
    
    [V,D] = eig(choi);
    D = max(real(D),0); % real is important, dtype complex spoils functioning of max. choi is Hermitian so D should be real. 
    choi = V*D*V';
    
    projected_choi_vec = reshape(choi,[],1);
end

