function [ projected_choi_vec ] = CPTNI_project( choi_vec)
%CPTNI_projectt: project a matrix into the set of CPTNI maps
% CPTNI stands for completely positive trace non-increasing
% choi_vec          : is a vector with dimensions (d^4 x 1).
%projected_choi_vec : (d^4 x 1) which represents a vectorised CPTNI Choi
%                   : matrix
    d = sqrt(size(choi_vec));
    
    x    = {choi_vec};
    GAP  = (1.0);
    k    = 1;
    p    = {0};
    q    = {0};
    y    = {0};
    while GAP(end) > 1e-12 
% %         k
%     for k=1:1000
%         x{k+1}=PSD_project(TNI_project(x{k}));     % ALTERNATING, seems to be fastest
%         x{k+1}=TNI_project(PSD_project(x{k}));     % alt - ALTERNATING, seems to be fastest
%         x{k+1}=0.5*PSD_project(x{k})+0.5*TNI_project(x{k});   % AVERAGED 
        y{k}   = TNI_project(x{k}+p{k}); % DIJKSTRA
        p{k+1} = x{k}+p{k}-y{k};
        x{k+1} = PSD_project(y{k}+q{k});
        q{k+1} = y{k}+q{k}-x{k+1};
        GAP    = norm(x{end-1}-x{end});  
        if k>2
            GAP(k)  = norm(p{k-1}-p{k})^2+norm(q{k-1}-q{k})^2+2*p{k-1}'*(x{k}-x{k-1}+2*q{k-1}'*(y{k}-y{k-1}));
%             GAP(k) = norm(q{k-1}-q{k})^2+norm(p{k-1}-p{k})^2;
%               GAP(k) = norm(x{k+1}-x{k});
        end
        k = k + 1;
   end
    projected_choi_vec = x{end};

end

