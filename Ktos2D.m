function [Rs] = Ktos2D(ncell,nei,ind,edges,cc,mid,flag)

% reconstruction operator from cells to edges

Rs = sparse(nei,ncell);

for e=1:nei
    K = edges(ind.internal(e),3);
    L = edges(ind.internal(e),4);
    if flag==1
        % centered
        Rs(e,K) = Rs(e,K) + 0.5;
        Rs(e,L) = Rs(e,L) + 0.5;
    elseif flag==2
        % linear reconstruction
        dK = sqrt((mid(ind.internal(e),1)-cc(K,1))^2+(mid(ind.internal(e),2)-cc(K,2))^2);
        dL = edges(ind.internal(e),5)-dK;
        Rs(e,K) = Rs(e,K) + dL/edges(ind.internal(e),5);
        Rs(e,L) = Rs(e,L) + dK/edges(ind.internal(e),5);
    else
        % mass average
        dK = sqrt((mid(ind.internal(e),1)-cc(K,1))^2+(mid(ind.internal(e),2)-cc(K,2))^2);
        dL = edges(ind.internal(e),5)-dK;
        Rs(e,K) = Rs(e,K) + dK/edges(ind.internal(e),5);
        Rs(e,L) = Rs(e,L) + dL/edges(ind.internal(e),5);
    end
end