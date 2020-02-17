function [JF] = JFk2D(ind,cc,Mx,grad,div,Rs,RK,RH,rhot,uk,tau,ddE)


ntri = size(cc,1);
nei = length(ind.internal);
phik = uk(1:ntri);
rhok = uk(ntri+1:2*ntri);
sk = uk(2*ntri+1:end);
ddEk = ddE(rhok,cc(:,1),cc(:,2));

gradphik = grad*phik;

% derivatives of the optimality conditions in phi
rho_all = [rhok;rhot];
I_all = [speye(ntri,ntri);zeros(ntri,ntri)];
JFpp = div*spdiags((Rs*RH*rho_all),0,nei,nei)*grad;
JFpr = Mx/tau+div*spdiags(gradphik,0,nei,nei)*Rs*RH*I_all;
JFps = sparse(ntri,ntri);

% derivatives of the optimality conditions in rho
JFrp = Mx/tau-(Mx*RH*I_all)'*RK*(spdiags(gradphik,0,nei,nei)*grad);
JFrr = Mx*spdiags(ddEk,0,ntri,ntri)/tau;
JFrs = -Mx/tau;

JFsp = sparse(ntri,ntri);
JFsr = Mx*spdiags(sk,0,ntri,ntri);
JFss = Mx*spdiags(rhok,0,ntri,ntri);

JF.pp = JFpp; JF.pr = JFpr; JF.ps = JFps;
JF.rp = JFrp; JF.rr = JFrr; JF.rs = JFrs;
JF.sp = JFsp; JF.sr = JFsr; JF.ss = JFss;

  

end