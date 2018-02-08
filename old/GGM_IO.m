function [ A ] = GGM_IO( d )
%GGM_IO preparations are Gellman eigenstates, measurements are Gellman
%expectations.
% the problem with this approach is that expectations can be negative, i.e.
% they are not probabilities! 
%   Detailed explanation goes here
    
% generate GGM inputs and outputs
    preparations = cell(1,d*(d*d-1));
    measurements = cell(1,d*d-1);
    i = 1; % indexes preparations
    m = 1; % indexes measurements
    for j=0:d-1
        for k=0:d-1
            if j==0 && k==0
                continue % don't need identity operator?? maybe do for trace decreasing
            end
            G = GenGellMann(j,k,d);
            [V,~] = eig(G);
            for l=1:d
                preparations{i} = V(:,l)*V(:,l)';
                i = i + 1;
            end
            measurements{m}=G;
            m = m + 1;
        end
    end
    
%     measurements = preparations;
    
    % construct GGM I/O matrix
    i = 1;
    num_measurements = length(measurements);
    num_preparations = length(preparations);
    A = zeros(num_measurements*num_preparations,d*d*d*d);
    for e=1:num_measurements
        for r=1:num_preparations
            E   = measurements{e};
            rho = preparations{r};
            row = reshape(kron(E',conj(rho)),[],1)'; % this was wrong before!
            A(i,:) = row;
            i = i+1;
        end
    end

    A = sparse(A);
    
end

