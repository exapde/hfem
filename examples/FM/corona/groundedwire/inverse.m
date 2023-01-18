function [rhoparam,UDG,UH] = inverse(mesh,master,app,uinf)

%(pi^(1/2)*c1*(erf(c2*(1/c3^2)^(1/2)) + erf((1/c3^2)^(1/2)*(2*pi - c2))))/(2*(1/c3^2)^(1/2))
%(sqrt(pi)/2)*(c1*c3)*(erf(c2/c3)+ erf((2*pi-c2)/c3)) 
%(uinf)/(2*pi*a*K*Ec)

% physical parameters
K = 1.5e-4; % (m^2/Vs)
e0 = 8.854e-12; % (F/m)

% lengthscale parameters  
a = 0.01; % radius of the wire  (m)
L0 = 1;   % reference length scale (m)  

Ec = 3.1*(1 + 0.0308*sqrt(1/a))*1e6; % Peek field
delta_peek = 0.15/100; % uncertainty in Peek field
Etmax = 40*1e3; % maximum background electric (V/m)
E0 = Etmax; % reference electric field (V/m)

% nondimensional parameters
Ecstar = Ec/E0;
Etstar = Etmax/E0;

% Estimate current per unit length based on the high wind approximation
% (Note: valid only for high wind)
Iwind = uinf*e0*E0;

% Estimate current per unit length when the ion density is set to one on the wire
Idens = 2*pi*a*K*(Ec+E0);

%
% 2*pi*a*K*(ub-lb)*1e6*(ub-lb)*1e6

% Estimate of the constant density (Note: valid only for high wind)
rhoe = Iwind/Idens;

% Normalization 
rho0 = rhoe*(L0/(e0*E0));

% Assume constant profile
rhoparam = [rho0 pi inf]
app.arg{10} = rhoparam;

% wind velocity field
[u,v] = potentialvelocity(mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:),uinf,a);
mesh.dgnodes(:,3,:) = u/(K*E0);
mesh.dgnodes(:,4,:) = v/(K*E0);

% Initial solution
UDG = initu(mesh,{1;rhoparam(1);0;0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver 
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Etstar; E = sqrt(Ex.^2+Ey.^2);
[dgnodes,u] = potentialfield(master,mesh,E,-5,1);    
[~,rho] = potentialfield(master,mesh,UDG(:,2,:),-5,1);    
x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
f = u(:,1,:);f=f(:);
q = rho(:,1,:);q=q(:);
t = cart2pol(x,y);
i = find(t<0);
t(i) = 2*pi+t(i);
[tj,jj] = sort(t);
[~,ii] = sort(abs(f-(1-delta_peek)*Ecstar));
[~,i1]=sort(t(ii(1:8)));
ii = ii(i1([1 end]));

% Determine the constant c2
[~,imax]=max(f(jj));
c2 = tj(imax);

% determine the constant c3
t1 = t(ii(1));
t2 = t(ii(2));
c3 =  min(abs(c2-t1),abs(t2-c2))/4;
 
% Estimate current per unit length when the ion density is a Gaussian
% distribution Exp(-(rho-c2)^2/c3^2)
Idens = a*K*(Ec+E0)*(sqrt(pi)/2)*(c3)*(erf(c2/c3)+ erf((2*pi-c2)/c3)); 
% Determine the constant c1
rhoe = Iwind/Idens;
c1 = rhoe*(L0/(e0*E0));

% Assume Gaussian profile c1*Exp(-(rho-c2)^2/c3^2)
rhoparam = [c1 c2 c3]
app.arg{10} = rhoparam;

% Initial solution
UDG = initu(mesh,{1;rhoparam(1);0;0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver 
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

Et = Etstar;     
Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
[dgnodes,u] = potentialfield(master,mesh,E,-5,1);    
[~,rho] = potentialfield(master,mesh,UDG(:,2,:),-5,1);    

x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
f = u(:,1,:);f=f(:);
q = rho(:,1,:);q=q(:);
t = cart2pol(x,y);
i = find(t<0);
t(i) = 2*pi+t(i);
[tj,jj] = sort(t);
[~,ii] = sort(abs(f-(1-delta_peek)*Ecstar));
[~,i1]=sort(t(ii(1:8)));
ii = ii(i1([1 end]));

figure(1);clf;
fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(f) min(f) max(f) max(f)]*E0/1e6,[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,E0*f(jj)/1e6,'b-',tj,E0*Ecstar*ones(size(t))/1e6,'--k','LineWidth',1.5);
plot(tj,(1+delta_peek)*E0*Ecstar*ones(size(t))/1e6,'-k','LineWidth',1.5);
plot(tj,(1-delta_peek)*E0*Ecstar*ones(size(t))/1e6,'-k','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));

figure(2);clf;
fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(q) min(q) max(q) max(q)]*e0*E0/L0*1e6,[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,1e6*e0*E0/L0*q(jj),'-b','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Charge density (\muC/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));
