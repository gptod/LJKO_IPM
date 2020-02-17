function grad = Grad2D(nei,ntri,ind,edges)


grad = sparse(nei,ntri);

for e=1:nei
    K = edges(ind.internal(e),3);
    L = edges(ind.internal(e),4);
    grad(e,K) = -1/edges(ind.internal(e),5);
    grad(e,L) = 1/edges(ind.internal(e),5);
end