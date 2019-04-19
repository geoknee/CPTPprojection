function [ projected_choi_vec ] = CPTP_project( choi_vec, MdagM, Mdagb  )
%CPTP_projectt: project a matrix into the set of CPTP maps
% CPTP stands for completely positive trace preserving
% choi_vec          : is a vector with dimensions (d^4 x 1).
% MdagM
% Mdagb             : helper matrices, previously computed 
%                   : which are used in the TP projection.
%projected_choi_vec : (d^4 x 1) which represents a vectorised CPTP Choi
%                   : matrix

    x    = {choi_vec};
    GAP  = (1);
    k    = 1;
    p    = {0};
    q    = {0};
    y    = {0};
    while GAP(end) >= 1e-4
        y{k}   = PSD_project(x{k}+p{k}); % DIJKSTRA
        p{k+1} = x{k}+p{k}-y{k};
        x{k+1} = TP_project(y{k}+q{k}, MdagM, Mdagb);
        q{k+1} = y{k}+q{k}-x{k+1}; 
        if k>2
            GAP(k)  = norm(p{k-1}-p{k})^2+norm(q{k-1}-q{k})^2+abs(2*p{k-1}'*(x{k}-x{k-1}))+abs(2*q{k-1}'*(y{k}-y{k-1})); 
        end
        k = k + 1;
   end
    projected_choi_vec = x{end};

end

