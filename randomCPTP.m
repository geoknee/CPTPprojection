function [ choi ] = randomCPTP(d, M)
%generate random completely positive, trace preserving maps in the Choi
%representation. 
%   following the first algorithm given in W. Bruzda and V. Cappellini and H.-J. Sommers and K. ?yczkowski
 %   Physics Letters A� 373� 320 - 324� (2009)
 % 


 X = randn(d*d,M)+ 1.0j*randn(d*d,M); %Ginibre matrix
 
 rho = X*X';
 
 Y = partial_trace(rho);
 
 LI = kron(inv(sqrtm(Y)),eye(d)); % our convention is different to Bruzda's
 
 choi = LI*rho*LI;
 
end

