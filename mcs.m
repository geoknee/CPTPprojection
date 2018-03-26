function [ cs ] = mcs(p,n )
%modified chi square
    
    cs = sum(((p-n).^2)./(p.*(1-p)));
    cs = sum(n)/length(n) * cs;
end

