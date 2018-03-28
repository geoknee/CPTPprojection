function [ cs ] = mcs(p,n,N )
%modified chi square
    
    cs = sum(((p-n).^2)./(p.*(1-p)));
    cs = 2* cs;
end

