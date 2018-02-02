function [ projected_choi_vec ] = PSD_project( choi_vec )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    d = sqrt(size(choi_vec));
    d = d(1);
    choi = reshape(choi_vec,[],d);
       
    choi = (choi + choi')/2;
    
    
    [V,D] = eig(choi);
    D = max(D,0);
    choi = V*D*V';

%     t = trace(choi);
%     choi = t*positiveProjection(choi/t); % Smolin method from EB. wants a
    %trace one matrix as input
    
    projected_choi_vec = reshape(choi,[],1);
end

