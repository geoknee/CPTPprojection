function [ A ] = PM_minimal( d )
%PM_minimal constructs a prepare/measure matrix featuring a minimal set of
% d^2 pure states
% and
% d^2 pairs of positive operators that sum to the identity 
% (i.e. a POVM with d^2 elements)
%. 2d^4 operators in all, each of size dxd.

    preparations = cell(1,d*d);
    i = 1;
    for k=1:d
        E = zeros(d,1);
        E(k)=1;
        preparations{i} = E*E';
        measurements{i} = E*E'/d^2;
        i = i+1;
    end
    
    for k=1:d-1
        for n=k+1:d
            E = zeros(d,1);
            E(k) = 1/sqrt(2);
            E(n) = 1/sqrt(2);
            preparations{i} = E*E';
            measurements{i} = E*E'/d^2;
            i = i+1;
        end
    end
    
    for k=1:d-1
        for n=k+1:d
            E = zeros(d,1);
            E(k) = 1/sqrt(2);
            E(n) = 1.0j/sqrt(2);
            preparations{i} = E*E';
            measurements{i} = E*E'/d^2; % to ensure proper normalisation
            i = i+1;
        end
    end
    
        
%     % make into a POVM by using complementary effects (each pair sums to
%     I/d^2 and there are d^2 such pairs)
    num_measurements = length(measurements);
    
    for i=1:num_measurements
        measurements{num_measurements+i} = eye(d)/d^2-measurements{i};
    end

    
    % construct I/O matrix
    i = 1;
    num_measurements = length(measurements);
    num_preparations = length(preparations);
    A = zeros(num_measurements*num_preparations,d*d*d*d);
    for r=1:num_preparations
        for e=1:num_measurements
            E   = measurements{e};
            rho = preparations{r};
            row = reshape(kron(conj(rho),E'),[],1)'; 

            A(i,:) = row;
            i = i+1;
        end
    end

    A = sparse(A);
    
end

