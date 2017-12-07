function [ out_list ] = additional_paulis( in_list )
%take a list of measurement operators and construct all tensor product
%combinations with a list of pauli operators
%   Detailed explanation goes here
paulis{1} = [1 0;0 0];
paulis{2} = [0 0;0 1];
paulis{3} = 0.5*[1 1;1 1];
paulis{4} = 0.5*[1 -1;-1 1];
paulis{5} = [0.5 -0.5j ;0.5j 0.5];
paulis{6} = [0.5 0.5j ;-0.5j 0.5];

l = length(in_list);

for i=1:l
    for j=1:6
        out_list{j+6*(i-1)} = kron(in_list{i},paulis{j});
    end
end

end

