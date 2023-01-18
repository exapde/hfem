setapplicationpath('FM/corona');

% See the paper Neimarlija2009.pdf and Medlin1998b.pdf
 
porder = 4;

% physical parameters
K = 1.5e-4; % (m^2/Vs)
D = 0.1;    %  (m^2/s)
e0 = 8.854e-12; % (F/m)

% lengthscale parameters  
a = 0.01; % radius of the wire  (m)
b = 15;   % distance from the wire to the plate (m)
L0 = 1;   % reference length scale (m)  

T0 = 1;   % reference time (s)
T1 = 20;  % time to reach maximum background electric (s)  
T2 = 30;  % final time (s)

% discharge electric strength
Ec = 3.1*(1 + 0.0308*sqrt(1/a))*1e6;
Etmax = 40*1e3; % maximum background electric (V/m)
E0 = Etmax; % reference electric field (V/m)

% nondimensional parameters
Kstar = K*E0*T0/L0;
Dstar = D*T0/L0;
estar = e0/e0;
Ecstar = Ec/E0;
Etstar = Etmax/E0;
T1star = T1/T0;
T2star = T2/T0;

% applied potential at the wire
phia = Etmax*b;  % (V)
phiastar = phia/(L0*E0);

% applied ion density at the wire
rhoa = 1e-9; % (C/m^3)
rhoastar = L0*rhoa/(e0*E0);

% stabilization parameter
beta1 = -1/200;
beta2  = -1/200;
tau = 10;
param = {Kstar,Dstar,estar,Ecstar,Etstar,T1star,T2star,beta1,beta2,tau};
ui = [phiastar rhoastar];

hybrid = 'hdg';
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fboucorona';
app.fhat = 'fhat';
app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
app.arg = param;
app.bcm = [6;8;7;9;1];
app.bcs = [[0 0];[0 0];[0 0];[0 0];ui];
app.bcd = [];
app.bcv = []; 

app.denseblock = 0;
app.tdep = true;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.np = 2;
app.nd = 2;
app.nch  = 2;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 2;

app.appname = 'corona';
app.time = [];
app.dtfc = [];
app.alpha = [];
app.dtcoef = [0,1]';

mesh = mkmesh_circleinrect2(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,app);

UDG = initu(mesh,{0;0;0;0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

wvel = wind(mesh,windvel);
mesh.dgnodes(:,3,:) = wvel/(K*E0);
mesh.dgnodes(:,4,:) = 0;

UDG0 = UDG; UH0 = UH;

nstage = 2;
torder = 2;
dt = 0.01;
ntime = 1930;
rhoaarray = zeros(ntime,1);
tarray = zeros(ntime,1);
time = 10.81;
savetime = 11:1:30;
facrho = 1.1;
facmul = 2;
for itime = 1:ntime
    fprintf('Timestep :  %d,   Time:  %g\n', itime, time);
       
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG0,UH0,[],time,dt,nstage,torder);        
    
    Et = min(((time+dt)/T1star),1.0)*Etstar;     
    Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
    [~,u] = potentialfield(master,mesh,E,-5,1);    
    Em = mean(abs(u(:)));  
    
    [Em Ecstar]
    if abs(Em-Ecstar)>0.0001*Ecstar
        % find lower bound and upper bound for rhoa
        if Em>Ecstar % keep increasing rhoa until Em<Ecstar
            while abs(Em-Ecstar)>0.0001*Ecstar
                rhol = rhoa; % Em>Ecstar
                rhoa = facmul*rhoa;                                                              
                rhoastar = L0*rhoa/(e0*E0);
                ui = [phiastar rhoastar];
                app.bcs = [[0 0];[0 0];[0 0];[0 0];ui];
                [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG0,UH0,[],time,dt,nstage,torder);    
                Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
                [~,u] = potentialfield(master,mesh,E,-5,1);    
                Em = mean(abs(u(:)));  
                if Em<Ecstar
                    rhou = rhoa; % Em<Ecstar
                    break;
                end
            end
        elseif Em<Ecstar % keep decreasing rhoa until Em>Ecstar
            while abs(Em-Ecstar)>0.0001*Ecstar
                rhou = rhoa; % Em<Ecstar
                rhoa = rhoa/facmul;                                                              
                rhoastar = L0*rhoa/(e0*E0);
                ui = [phiastar rhoastar];
                app.bcs = [[0 0];[0 0];[0 0];[0 0];ui];
                [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG0,UH0,[],time,dt,nstage,torder);    
                Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
                [~,u] = potentialfield(master,mesh,E,-5,1);    
                Em = mean(abs(u(:)));  
                if Em>Ecstar
                    rhol = rhoa; % Em<Ecstar
                    break;
                end
            end
        end
        while abs(Em-Ecstar)>0.0001*Ecstar                
            rhoa = 0.5*(rhol+rhou);                    
            rhoastar = L0*rhoa/(e0*E0);
            ui = [phiastar rhoastar];
            app.bcs = [[0 0];[0 0];[0 0];[0 0];ui];
            [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG0,UH0,[],time,dt,nstage,torder);    
            Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
            [~,u] = potentialfield(master,mesh,E,-5,1);    
            Em = mean(abs(u(:)));  
            if Em<Ecstar % need to increase Em -> reduce rhoa
                rhou = rhoa;
            else % need to decrease Em -> increase rhoa
                rhol = rhoa;
            end
        end                    
    end
    UDG0 = UDG; UH0 = UH;        
                
    if time <T1star        
        rhoa=facrho*rhoa;
        rhoastar = L0*rhoa/(e0*E0);
        ui = [phiastar rhoastar];
        app.bcs = [[0 0];[0 0];[0 0];[0 0];ui];
    end        
    
    time = time + dt;    
    tarray(itime) = time;
    rhoaarray(itime) = rhoa;

    if itime>3
        facrho = rhoaarray(itime)/rhoaarray(itime-1);
        facmul = max(facrho,1.2);
    end
    
    [Em Ecstar]        
    rhoaarray(1:itime)'
    
    for j=1:length(savetime)        
        if abs(time-savetime(j))<1e-10
            fn = ['solc' num2str(savetime(j)) 'wind' num2str(windvel) '.mat'];
            save(fn,'UDG','UH','rhoaarray','tarray');
        end
    end    
end
figure(1); clf; scaplot(mesh,(L0*E0)*UDG(:,1,:),[],2); axis equal; axis tight; colormap jet; 
figure(2); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:),[],2); axis equal; axis tight; colormap jet;    


return;

% r = sqrt(mesh.dgnodes(:,1,:).^2 + mesh.dgnodes(:,2,:).^2);
% UDG(:,1,:) = b-r;
% UDG(:,2,:) = 1;
% UDG(:,3,:) = mesh.dgnodes(:,1,:)./r;
% UDG(:,5,:) = mesh.dgnodes(:,2,:)./r;
% UDG = exactsol(mesh.dgnodes, cell2mat(param));
% UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

for i = 1:2
    if i==1
        tm = L0*E0;
    else
        tm = e0*E0/L0;
    end
    figure(i); clf; 
    scaplot(mesh,tm*UDG(:,i,:),[],2); 
    axis equal; axis tight; axis on; colormap jet;
    colorbar('FontSize',15);
    set(gca,'FontSize',18);
    xlabel('x (m)','FontSize',20);
    ylabel('y (m)','FontSize',20);
    box off;
    %set(gca, 'LooseInset', get(gca, 'TightInset'));
    fn = ['numsol' num2str(i)];
    print('-dpng',fn);
end

Ex = E0*UDG(:,3,:);
Ey = E0*UDG(:,5,:);
E = sqrt(Ex.^2+Ey.^2);
figure(3); clf; 
scaplot(mesh,E,[],2); 
axis equal; axis tight; axis on; colormap jet;
colorbar('FontSize',15);
set(gca,'FontSize',18);
xlabel('x (m)','FontSize',20);
ylabel('y (m)','FontSize',20);
box off;

%set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['numsol' num2str(3)];
print('-dpng',fn);

return;

[dgnodes,u] = potentialfield(master,mesh,E,-1,[1]);
figure(4); clf;
hold on;
for k = 1:size(u,3)    
    plot(dgnodes(:,1,k),u(:,1,k),'-k','LineWidth',1.5);
end
hold off;
axis normal; axis tight; axis on; box on;
xlabel('Distance from centerline (m)','FontSize',20);
ylabel('Electric field intensity (V/m)','FontSize',20);
set(gca,'xtick',[-4:1:8]);
%set(gca,'ytick',[0:0.2:1]);
axis([-4 8 0 2e5])
set(gca,'FontSize',18);


J(:,1,:) = K*(e0*E0/L0)*E0*UDG(:,2,:).*(UDG(:,3,:) + mesh.dgnodes(:,3,:));
J(:,2,:) = K*(e0*E0/L0)*E0*UDG(:,2,:).*(UDG(:,5,:) + mesh.dgnodes(:,4,:));
Ja = sqrt(J(:,1,:).^2 + J(:,2,:).^2);
[dgnodes,u] = potentialfield(master,mesh,Ja,-1,[1]);
figure(5); clf;
hold on;
for k = 1:size(u,3)    
    plot(dgnodes(:,1,k),1e6*u(:,1,k),'-k','LineWidth',1.5);
end
hold off;
axis normal; axis tight; axis on; box on;
xlabel('Distance from centerline (m)','FontSize',20);
ylabel('Current density (\muA/m^2)','FontSize',20);
axis([-4 6 0 15])
set(gca,'xtick',[-4:1:6]);
set(gca,'ytick',[0:1:15]);
set(gca,'FontSize',18);


yc = -1.955412550450295;
%yc = -1.915589847848507;
x  = mesh.dgnodes(:,1,:);
y  = mesh.dgnodes(:,2,:);
in = find(abs(y+2)<1e-6); 
xc = x(in);
jc = Ja(in);
figure(7); clf; 
plot(xc,1e6*jc,'o');
axis normal; axis tight; axis on; box on;
xlabel('Distance from centerline (m)','FontSize',20);
ylabel('Current density (\muA/m^2)','FontSize',20);
axis([-4 8 0 15])
set(gca,'xtick',[-4:1:8]);
set(gca,'ytick',[0:1:15]);
set(gca,'FontSize',18);









% 
% meshb.p(:,1) = -mesha.p(:,1);
% figure(1); clf; 
% meshplot(mesha); 
% hold on;
% meshplot(meshb); 
% axis equal; axis tight; axis off; 
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% 





