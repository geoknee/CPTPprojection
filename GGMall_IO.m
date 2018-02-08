function [ A ] = GGMall_IO( d )
%GGMall_IO constructs a prepare/measure matrix featuring all eigenstates of
%all but the trivial GGM matrices. so total number of combinations of
%operators is d*(d*d-1) squared: d^6 overall.
    
% generate GGM inputs and outputs
    preparations = cell(1,d*(d*d-1));
    i = 1;
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
        end
    end
    
    measurements = preparations;
    
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

