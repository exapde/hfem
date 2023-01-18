setapplicationpath('FM/corona');

% See the paper Neimarlija2009.pdf and Medlin1998b.pdf

% fn = ['solfwind' num2str(windvel) '.mat']; 
% load(fn);

% Etmax, delta_peek, windlevel

porder = 6;

% physical parameters
K = 1.5e-4; % (m^2/Vs)
D = 0.0005;    %  (m^2/s)
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
%Etmax = 40*1e3; % maximum background electric (V/m)
E0 = Etmax; % reference electric field (V/m)
%delta_peek = 0.1/100;

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
%rhoa = 0.045775*1.5796e-04; % (C/m^3)
% rhoastar = L0*rhoa/(e0*E0);
% rhoparam = [L0*rhoa/(e0*E0) 0.9835*pi 0.2495];

lb =  E0*Ecstar/1e6;
ub = (1+delta_peek)*lb;
% if windvel==200
%     beta = E0*windvel/(K*a*pi*(ub^2-lb^2)*1e6^2);
% end
% beta200=4.128125367460415e+02;
% beta100=1.935058765997070e+02;
% beta50 =  87.957216636230442;
% if windvel>100
%     beta = beta200*(windvel/200);
% elseif windvel>50
%     beta = beta100*(windvel/100);
% else
%     beta = beta50*(windvel/50);
% end
%beta = E0*windvel/(0.25*K*a*pi*(ub^2-lb^2)*1e6^2);

% stabilization parameter
beta1 = -1/100;
beta2  = -1/100;
rhoparam = [beta Ecstar];
tau = max(10,windvel);
param = {Kstar,Dstar,estar,Ecstar,Etstar,T1star,T2star,beta1,beta2,rhoparam,tau};
ui = [phiastar rhoparam(1)];

hybrid = 'hdg';
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fboucorona2';
app.fhat = 'fhat';
app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
app.arg = param;
if windvel>0
    app.bcm = [10;8;8;9;0];
else
    app.bcm = [10;9;8;9;0];
end 
app.bcs = [[0 0];[0 0];[0 0];[0 0];ui];
app.bcd = [];
app.bcv = []; 

app.denseblock = 0;
app.tdep = false;
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

mesh = mkmesh_circleinrect5(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,app);

[u,v] = potentialvelocity(mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:),windvel,0.01);
mesh.dgnodes(:,3,:) = u/(K*E0);
mesh.dgnodes(:,4,:) = v/(K*E0);

% if windvel==200
%     UDG = initu(mesh,{1;rhoparam(1);0;0;0;0});
%     UH = inituhat(master,mesh.elcon,UDG,app.ncu);
% end

[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
plotsol;

[~,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
ewire = u(:,end,:)*E0/1e6;
maxe = max(abs(ewire(:)));    
[lb ub maxe]
    
% determine betalb and betaub
tol = 0.02;
if (abs(maxe-ub)/(ub-lb) < tol)        
    return;
elseif maxe>ub
    betalb = beta;        
    while (1)
        beta = 2*beta;
        app.arg{10}(1) = beta;    
        [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
        plotsol;
        [~,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
        ewire = u(:,end,:)*E0/1e6;
        maxe = max(abs(ewire(:)));      
        if (abs(maxe-ub)/(ub-lb) < tol)        
            break;
        elseif maxe>ub
            betalb = beta;            
        else
            betaub = beta;
            break;
        end
        [lb ub maxe]
    end
elseif maxe<ub
    betaub = beta;
    while (1)
        beta = beta/2;
        app.arg{10}(1) = beta;    
        [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
        plotsol;
        [~,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
        ewire = u(:,end,:)*E0/1e6;
        maxe = max(abs(ewire(:)));        
        if (abs(maxe-ub)/(ub-lb) < tol)        
            break;
        elseif maxe<ub
            betaub = beta;            
        else
            betalb = beta;
            break;
        end
        [lb ub maxe]
    end    
end

[~,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
ewire = u(:,end,:)*E0/1e6;
maxe = max(abs(ewire(:)));        
while (abs(maxe-ub)/(ub-lb) > tol)     
    [lb ub maxe]
    [betalb betaub]
    beta = 0.5*(betalb+betaub);
    app.arg{10}(1) = beta;    
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
    plotsol;
    [~,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
    ewire = u(:,end,:)*E0/1e6;
    maxe = max(abs(ewire(:)));      
    if maxe>ub
        betalb = beta;            
    else
        betaub = beta;    
    end                
end    


%     app.arg{10}(1) = app.arg{10}(1)*(1 + (maxe-ub)/(ub-lb));
%     app.arg{10}
%     [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
%     [~,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
%     ewire = u(:,end,:)*E0/1e6;
%     maxe = max(abs(ewire(:)));        
%     plotsol;

