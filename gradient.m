function [ g ] = gradient( A, n ,choi_vec)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    p   = A*choi_vec;
    p   = abs(p/sum(p));
    eta = n./p;
%     g   = -(eta.'*A).';
    g = -A'*eta;
end
