function [ d ] = trace_dist( A,B )
%trace_dist calculates the trace distance between density matrixes A and B
    d = 0.5*sum(svd(A-B));
end

