function [F] = Fk2D(ind,cc,Mx,grad,div,Rs,RK,RH,rhot,uk,tau,dE,mu)

ntri = size(cc,1);
nei = length(ind.internal);
phik = uk(1:ntri);
rhok = uk(ntri+1:2*ntri);
sk = uk(2*ntri+1:end);
dEk = dE(rhok,cc(:,1),cc(:,2));

gradphik = grad*phik;

% optimality conditions in phi
rho_all = [rhok;rhot];
Fp = Mx*(rhok-rhot)/tau+div*((Rs*RH*rho_all).*gradphik);

% optimality conditions in rho
I_all = [speye(ntri,ntri);zeros(ntri,ntri)];
Fr = Mx*(phik+dEk-sk)/tau-0.5*(Mx*RH*I_all)'*RK*(gradphik).^2;

Fs = Mx*(rhok.*sk-mu);

F.p = Fp;
F.r = Fr;
F.s = Fs;


end

