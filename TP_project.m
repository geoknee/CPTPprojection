    function [ projected_choi_vec ] = TP_project( choi_vec )
    %UNTITLED6 Summary of this function goes here
    %   Detailed explanation goes here
        d = sqrt(sqrt(size(choi_vec)));
        d = d(1);
        % initialise A. TODO precompute this!
        M = zeros([d*d,d*d*d*d]);
        for i=1:d
            e = zeros(1,d);
            e(i)  = 1;
            B = kron(eye(d),e);
            M = M + kron(B,B);
        end
        b = reshape(eye(d),[],1);
        S = eye(d*d*d*d)- (1.0/d) * (M'*M);
        projected_choi_vec = S*choi_vec + (1.0/d)*M'*b;
    end

