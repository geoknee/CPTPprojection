function [ A ] = PM_minimal( d )
%PM_minimal constructs a prepare/measure matrix featuring a minimal set of projectors
%. so total number of combinations of
% operators is d*d squared: d^4 overall.
%   Detailed explanation goes here
    
% generate GGM inputs and outputs
    preparations = cell(1,d*d);
    i = 1;
    for k=1:d
        E = zeros(d,1);
        E(k)=1;
        preparations{i} = E*E';
        i = i+1;
    end
    
    for k=1:d-1
        for n=k+1:d
            E = zeros(d,1);
            E(k) = 1/sqrt(2);
            E(n) = 1/sqrt(2);
            preparations{i} = E*E';
            i = i+1;
        end
    end
    
    for k=1:d-1
        for n=k+1:d
            E = zeros(d,1);
            E(k) = 1/sqrt(2);
            E(n) = 1.0j/sqrt(2);
            preparations{i} = E*E';
            i = i+1;
        end
    end
    
    
    measurements = preparations;
    
    % construct I/O matrix
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

