for d=2:10
    d
    A = PM_minimal(d);
    density_list(d) = nnz(A)/prod(size(A));
end
figure
plot(density_list(2:end))
xlabel('d')
ylabel('density of A')
saveas(gcf,'density.png')
saveas(gcf,'density.eps','epsc')
