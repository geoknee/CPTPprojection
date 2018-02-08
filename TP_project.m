function [ projected_choi_vec ] = TP_project( choi_vec, MdagM, Mdagb)
%TP_projectt: project a matrix into the set of trace preserving maps
% choi_vec          : is a vector with dimensions (d^4 x 1).
% MdagM
% Mdagb             : helper matrices, previously computed, 
%                   : which are used in the TP projection.
%projected_choi_vec : (d^4 x 1) which represents a vectorised TP Choi
%                   : matrix
    d = sqrt(sqrt(size(choi_vec)));
    d = d(1);
    
%     
%     S = sparse(eye(d*d*d*d)- (1.0/d) * MdagM);
%     projected_choi_vec = S*choi_vec + (1.0/d)*Mdagb;
    
    projected_choi_vec = choi_vec -(1/d)*MdagM*choi_vec + (1.0/d)*Mdagb;
    
end
