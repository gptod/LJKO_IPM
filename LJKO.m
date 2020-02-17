clear
close all

% load mesh
load('sq_mesh5')

%% DEFINITION OF THE GRADIENT FLOW:

% Fokker-Planck equation:
% definition of the energy
% g = 1;
% V =@(x,y) -g*x;
% E =@(rho,x,y) rho.*log(rho)+rho.*V(x,y);
% dE =@(rho,x,y) log(rho)+1+V(x,y);
% ddE =@(rho,x,y) 1./rho;
% % definition of the exact solution
% lx = 1;
% g = 1;
% a = lx*(pi^2+0.25*g^2);
% u =@(x,y,t) exp(-a*t+0.5*g*x).*(pi*cos(pi*x)+0.5*g*sin(pi*x))+pi*exp(g*(x-0.5));

% Porous media equation:
% definition of the energy
V =@(x,y) 0.5*((x-0.5).^2+(y-0.5).^2);
m = 2; % m > 1
E =@(rho,x,y) rho.^m/(m-1)+rho.*V(x,y); 
dE =@(rho,x,y) m*rho.^(m-1)/(m-1)+V(x,y);
ddE =@(rho,x,y) m*rho.^(m-2);

% Barenblatt profile for the equilibrium:
% BB =@(x,y,mass) max((0.5*mass/pi)^((m-1)/m)-((m-1)/(2*m))*((x-0.5).^2+(y-0.5).^2),0).^(1/(m-1));
% rhoBB = BB(cc(:,1),cc(:,2),mass);
% EBB = sum(area.*E(rhoBB,cc(:,1),cc(:,2)));


%% INITIAL CONDITIONS:

% exact solution for the F-P
% rho0 = u(cc(:,1),cc(:,2),0);

% gaussian density
initial =@(x,y) exp((-(x-0.5).^2-(y-0.5).^2)/0.01);
rho0 = initial(cc(:,1),cc(:,2));

% % delta density
% initial =@(x,y) 100*(0.47<=x).*(x<=0.53).*(0.47<=y).*(y<=0.53);
% rho0 = initial(cc(:,1),cc(:,2));

% % cross shaped initial condition
% rot = [sqrt(2)/2  sqrt(2)/2; - sqrt(2)/2  sqrt(2)/2];
% initial = @(x,y) (abs(x-0.5)<0.1).*(abs(y-0.5)<0.4)+(abs(y-0.5)<0.1).*(abs(x-0.5)<0.4);
% prot = ([cc(:,1) cc(:,2)]-0.5)*rot+0.5;
% rho0 = initial(prot(:,1),prot(:,2)); rho0(rho0>0) = rho0(rho0>0)./rho0(rho0>0);

% total mass
mass = sum(area.*rho0);
%rho0 = rho0./mass; mass=sum(area.*rho0); % normalization


%% PARAMETERS

% choose reconstruction: 1 -> 0.5, 2 -> linear, 3 -> mass                      
flag = 3;

T = 1; % integration time
tau0 = 0.05; % initial time step
tau_max = 0.1; % maximum time step
eps_0 = 1e-6; % tolerance

eloc = 1; % set 1 if the energy is local, 0 if not
verb = 2; % verbosity level: {0,1,2}

% IPM scheme's parameters:
k1max = 20; %minimum number of outer iteration
k2max = 10; %maximum number of inner iteration
gamma0 = 0.1; % decay ratio


%% ASSEMBLE MATRICES
Mx = spdiags(area,0,ncell,ncell); % mass matrix for cells
Ms = spdiags(edges(ind.internal,5).*edges(ind.internal,6),0,nei,nei); % mass matrix for internal edges
RH = assembleRH(ncell); % weight in time matrix
grad = Grad2D(nei,ncell,ind,edges); % gradient matrix
div = -(Ms*grad)'; % divergence matrix
Rs = Ktos2D(ncell,nei,ind,edges,cc,mid,flag); % reconstruction operator on the edges
RK = Mx\(Rs'*Ms); % reconstruction operator on the cells


%% INITIALIZATION

phit = dE(rho0,cc(:,1),cc(:,2));
rhot = rho0;
rhotn = max(5e-2,rho0); rhotn = rhotn*mass/sum(area.*rhotn);
%rhotn = mass./area;
stn = ones(ncell,1);
uk = [phit; rhotn; stn];

time = 0;
itn = 0;

Energy = sum(E(rho0,cc(:,1),cc(:,2)).*area);
ts = time;


%% TIME INTEGRATION

while abs(time-T) > 1e-10
    
    tau = tau0;
    time = time+tau;
    if time>T && abs(time-T)>1e-10
        time = time-tau;
        tau = T-time;
        time = time+tau;
    end
    itn = itn+1;

    itk1 = 0;
    itk2 = 0;
    errs = [];
    tit = 0;
    mu0 = sum(uk(ncell+1:2*ncell).*uk(2*ncell+1:3*ncell))/ncell;
    gamma = gamma0;
    mu = mu0/gamma;
    
    % outer cycle
    while true
        
        Fk = Fk2D(ind,cc,Mx,grad,div,Rs,RK,RH,rhot,uk,tau,dE,0);
        delta_0 = norm([Fk.p;Fk.r;Fk.s],2);
        errs = [errs; delta_0];
        if verb>0
            fprintf('%12s %4i %7s %1.4e %19s %4i %10s %1.4e \n','Outer step: ',itk1,'Error: ',delta_0,'Newton iterations: ',itk2,'rho min: ',min(uk(ncell+1:2*ncell)))
        end
        if delta_0 < eps_0
            break
        end
        itk1 = itk1+1;
        if itk1 > k1max
            time = time-tau;
            tau0 = tau0/2;
            tau = tau0;
            time = time+tau;
            uk(1:ncell) = dE(rhot,cc(:,1),cc(:,2));
            uk(ncell+1:2*ncell) = rhotn;
            uk(2*ncell+1:end) = ones(ncell,1);
            itk1 = 0;
            itk2 = 0;
            errs = [];
            mu0 = sum(uk(ncell+1:2*ncell).*uk(2*ncell+1:3*ncell))/ncell;
            gamma = gamma0;
            mu = mu0/gamma;
            fprintf('%10s %4i %24s %1.5e \n','iteration ',itn,', decrease time step to',tau)
            continue
        end
        mu = gamma*mu;
        eps_mu = eps_0;
        
        %inner cycle
        itk2 = 0;
        while true
            
            Fk = Fk2D(ind,cc,Mx,grad,div,Rs,RK,RH,rhot,uk,tau,dE,mu);
            delta_mu = norm([Fk.p;Fk.r;Fk.s],2);
            if verb>1
               fprintf('%12s %4i %8s %1.4e \n','Inner step: ',itk2,'Error: ',delta_mu)
            end
            if delta_mu < eps_mu
                if itk1==1
                    if itk2<8 && 2*tau0<=tau_max
                        tau0 = 2*tau0;
                    end
                end
                tit = tit+itk2;
                break
            end
            itk2 = itk2+1;
            if itk2 > k2max
                if itk1==1
                    itk1=k1max;
                    break
                end
                mu = mu/gamma;
                gamma = gamma+0.3*(1-gamma);
                mu = gamma*mu;
                uk = [phimu;rhomu;smu];
                tit = tit+itk2;
                itk2 = 0;
                continue 
            end
            
            JFk = JFk2D(ind,cc,Mx,grad,div,Rs,RK,RH,rhot,uk,tau,ddE);
            A = JFk.pp; B = JFk.pr; C = JFk.rp; D = JFk.rr-JFk.rs*(JFk.ss\JFk.sr);
            f1 = Fk.p;
            f2 = Fk.r-JFk.rs*(JFk.ss\Fk.s);
            if eloc==1
                S = A-B*(D\C);
                dp = -S\(f1-B*(D\f2));
                dr = D\(-f2-C*dp);
                ds = JFk.ss\(-Fk.s-JFk.sr*dr);
            else
                d = [A B; C D]\[-f1;-f2]; dp = d(1:ncell); dr = d(ncell+1:2*ncell);
                ds = JFk.ss\(-Fk.s-JFk.sr*dr);
            end
            
            alfak = 1;
            while any(uk(ncell+1:2*ncell)+alfak*dr<=0)
                alfak = 0.5*alfak;
            end
            uk = uk + alfak*[dp;dr;ds];
                        
        end
        
        phimu = uk(1:ncell);
        rhomu = uk(ncell+1:2*ncell);
        smu = uk(2*ncell+1:end);
        
    end
    
    fprintf('%11s %4i %22s %4i %7s %1.4e \n','Time step: ',itn,'Total Newton iterations: ',tit,'Error: ',delta_0)
    
    phit = uk(1:ncell);
    rhot = uk(ncell+1:2*ncell);
    rhotn = max(1e-2,rhot); rhotn=rhotn*mass/sum(area.*rhotn);
    uk(ncell+1:2*ncell) = rhotn;
    uk(2*ncell+1:end) = ones(ncell,1);
    
    Energy = [Energy; sum(E(rhot,cc(:,1),cc(:,2)).*area)];
    ts = [ts; time];
    
    Frho = scatteredInterpolant(cc,rhot);
    Zrho = Frho(nodes(:,1),nodes(:,2));
    figure(3)
    trisurf(cells(:,2:end),nodes(:,1),nodes(:,2),Zrho)
    %axis([0 1 0 1 0 1.2])
    colormap('jet')
  

end




