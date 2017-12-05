choi = rand(4)+1.0j*rand(4);
choi = (choi + choi')/2;
eig(choi)
choi
for i=1:1000
    choi = (choi + choi')/2;
    [V,D] = eig(choi);
    D = max(D,0);
    choi = V*D*V';
end
eig(choi)
choi