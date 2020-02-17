function RH = assembleRH(ncell)

% recontruction of the weight density for the H^-1 norm

% final density as weight in the H^{-1} norm
RH = [speye(ncell,ncell) sparse(ncell,ncell)];
% initial density as weight in the H^{-1} norm
%RH = [sparse(ntri,ntri) speye(ntri,ntri)];
% arithmetic mean as weight in the H^{-1} norm
%RH = [0.5*speye(ntri,ntri) 0.5*speye(ntri,ntri)];
