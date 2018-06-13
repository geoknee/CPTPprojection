function [ d ] = trace_dist( A,B )
%trace_dist calculates the trace distance between density matrixes A and B

    if isnan(A-B)
        d = NaN;
    else
        d = 0.5*sum(svd(A-B));
    end
    
    
    %debugging, use frobenius-norm instead
%     d = norm(A-B,'fro');
%     d = 0;
end

