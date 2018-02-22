function [ projected_choi_vec ] = PSD_project( choi_vec )
%PSD_projectt: project a matrix into the set of completely positive maps
% we rely on Choi's theorem, that the Choi matrix's poitivty <--> the
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
    D = max(D,0);
    choi = V*D*V';

%     t = trace(choi);
%     choi = t*positiveProjection(choi/t); % Smolin method from EB. wants a
    %trace one matrix as input
    
    projected_choi_vec = reshape(choi,[],1);
end

