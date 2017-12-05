function [ c ] = cost( A, n , choi_vec )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    p = A*choi_vec;
    p = abs(p/sum(p));
    c = real(-(n.'*reallog(p)));
end

