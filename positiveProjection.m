function [rhoOut]= positiveProjection(rho)

d=length(rho);

% rho=rho/trace(rho);

[V,Dunsort] = eig(rho);

Dun = real(diag(Dunsort));
[D,I] = sort(Dun,1,'descend');
V = V(:, I);

lambda = zeros(d,1);
a=0;
for i=d:-1:1
    if D(i) + a/i >= 0
        for j=1:i
            lambda(j) = D(j) + a/i;
        end
        break
    else
        lambda(i) = 0;
        a = a + D(i);
    end
end

Dl = diag(lambda);
rhoOut = V*Dl*V';



