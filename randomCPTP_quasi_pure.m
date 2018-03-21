function [ choi ] = randomCPTP_quasi_pure(d, purity)
%generate random completely positive, trace preserving maps in the Choi
%representation with full rank and purity p
%  


lambda=0;
purityTemp=0;
%Generate exponialy decreasing eigenvalues with the specified purity
while purityTemp<purity
    lambda = lambda + 0.001; %increase std until reach correct purity
    lam = real(exp(-lambda*(1:d*d))); %exponential distribution of eigenvalues
    lamb = lam/sum(lam);
    purityTemp = sum(lamb.^2);
end

choi = 0;
for i=1:d*d
    choi = choi + lamb(i)*randomCPTP(d,1);
end

choi = d*choi/trace(choi);

% trace(choi*choi)/(d*d);

