function [ d ] = trace_dist( A,B )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [~,S,~] = svd(A-B);
    d = 0.5*trace(S);
    % nicely indented
end

