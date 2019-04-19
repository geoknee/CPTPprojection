function [ rho_reduced ] = partial_trace( rho )
%rho_reduced computes the partical trace of a density matrix
%   Detailed explanation goes here
  d = sqrt(size(rho));
  d = d(1);
  rho_reduced = 0;
    for i=1:d
        e = zeros([1,d]);
        e(i) = 1;
        bra = kron(eye(d),e);
        ket = bra';
        rho_reduced = rho_reduced + bra*rho*ket;
    end
end

