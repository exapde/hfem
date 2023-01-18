function [UDGEX,UDGIM,UHEX,UHIM,PDGEX,PDGIM] = hdg_solve_IMEX(master,meshEX,meshIM,app,UDGEX,UDGIM,UHEX,UHIM,PDGEX,PDGIM,Minv,time,dt,nstage)
%HDG_SOLVE_IMEX Solve using the HDG method and Implicit-Explicit DIRK timestepping.
% 
% SYNTAX: [UDGEX,UDGIM,UHEX,UHIM] = HDG_SOLVE_IMEX(MASTER,MESHEX,MESHIM,APP,UDGEX,UDGIM,UHEX,UHIM,MINV,TIME,DT,NSTAGE)
%
% INPUTS:
%    MASTER:                  Master structure
%    MESHEX:                  Mesh structure for Explicit domain
%    MESHIM:                  Mesh structure for Implicit domain
%    APP:                     Application structure
%    UDGEX(npv,nc,neEX):      Vector of unknowns, explicit domain (initial guess)
%    UDGIM(npv,nc,neIM):      Vector of unknowns, implicit domain (initial guess)
%    UHEX(nch,npf*nfEX):      Vector of hybrid unknowns, explicit domain (initial guess)
%    UHIM(nch,npf*nfIM):      Vector of hybrid unknowns, implicit domain (initial guess)
%    PDGEX(npv,nc,neEX):      Vector of time integrated variables, explicit domain (initial guess)
%    PDGIM(npv,nc,neIM):      Vector of time integrated variables, implicit domain (initial guess)
%    MINV:                    Inverted mass matrix
%    TIME:                    Time
%    DT:                      Timestep
%    NSTAGE:                  Number of stages
%
% OUTPUTS:
%    UDGEX(npv,nc,neEX):      Vector of unknowns, explicit domain (updated)
%    UDGIM(npv,nc,neIM):      Vector of unknowns, implicit domain (updated)
%    UHEX(nch,npf*nfEX):      Vector of hybrid unknowns, explicit domain (updated)
%    UHIM(nch,npf*nfIM):      Vector of hybrid unknowns, implicit domain (updated)
%    PDGEX(npv,nc,neEX):      Vector of time integrated variables, explicit domain (updated)
%    PDGIM(npv,nc,neIM):      Vector of time integrated variables, implicit domain (updated)
%
% VECTORS DIMENSIONS
%    npv:                     Number of DG nodes within an element
%    nc:                      Number of conservation equations solved (components)
%    nch:                     Number of components for the hybrid variable
%    npf:                     Number of HDG nodes per face
%    neEX:                    Number of elements in the explicit domain
%    neIM:                    Number of elements in the implicit domain
%    nfEX:                    Number of faces in the explicit domain
%    nfIM:                    Number of faces in the implicit domain
%
% Ref: High-Order LES Simulations using Implicit-Explicit Runge-Kutta 
% Schemes Diagonally Implicit Runge Kutta Methods for Stiff ODE's, 
% by Per-Olof Persson, 49th AIAA Aerospace Sciences Meeting, January 2011.
% 
% Author(s): Sebastien Terrana, Lauren Kolkman
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% March 2018

% Global IMEX data
global IMEXdata;

% Get DIRK coefficients
[Ainv,bAinv,c] = dirkcoeff(nstage,nstage);
Ainv = Ainv/dt;

% Get ERK coefficients
[Ahat,bhat,chat] = ERKcoeff(nstage+1,nstage);
% Values at begining of time step
UDGEXn = UDGEX;
UHEXn  = UHEX;
PDGEXn = PDGEX;
UDGIMn = UDGIM;
PDGIMn = PDGIM;
UHIMn  = UHIM;

% Allocate Residuals
[npv,ncu,ne] = size(UDGEX);
ResEX = zeros(npv,ncu,ne,nstage+1);
ResEXpdg = zeros(npv,1,ne,nstage+1);
% ResIM = zeros(npv,ncu,ne,nstage);
UDGIM = (1-sum(bAinv))*UDGIMn;
PDGIM = (1-sum(bAinv))*PDGIMn;
UHIM  = (1-sum(bAinv))*UHIMn;

% Transfer data to Explicit domain
IMEXdata = getVhatIMEX(meshIM,UHIMn);
% First Explicit Residual
ResEX(:,:,:,1) = hdgres(master,meshEX,app,UDGEX);
ResEXpdg(:,:,:,1) = UDGEXn(:,1,:);

for istage = 1:nstage
    
    % Assign Diagonal DIRK coeff
    app.fc_u = Ainv(istage,istage);
    app.fc_q = Ainv(istage,istage);
    
    % Build current stage second member
    ResSum = 0;
    ResSumPdg = 0;
    for j=1:istage
        ResSum = ResSum + Ahat(istage+1,j)*ResEX(:,:,:,j);
        ResSumPdg = ResSumPdg + Ahat(istage+1,j)*ResEXpdg(:,:,:,j);
    end
    
    % Update Explicit Solution : UDGEX = UDGEXn - dt*Minv*ResSum;
    for i = 1:ne
        UDGEX(:,:,i) = UDGEXn(:,:,i) - dt*Minv(:,:,i)*ResSum(:,:,i); 
        PDGEX(:,:,i) = PDGEXn(:,:,i) + dt*ResSumPdg(:,:,i);   
    end
    
    % Transfer data to Implicit domain
    IMEXdata = getQhatIMEX(master,meshEX,app,UDGEX);
    
    % Source Term for the implicit solve
    if istage==1
        SDG = Ainv(1,1)*UDGIMn; 
    elseif istage == 2
        SDG = (Ainv(2,1)+Ainv(2,2))*UDGIMn - Ainv(2,1)*UDGIMtmp;
    elseif istage == 3
        SDG = (Ainv(3,1)+Ainv(3,2)+Ainv(3,3))*UDGIMn - ...
              (Ainv(3,1)/Ainv(2,1))*((Ainv(2,1)+Ainv(2,2))*UDGIMn-SDG) - Ainv(3,2)*UDGIMtmp;    
    end
    
    % Update local time
    app.time = time+dt*c(istage);
    
    % Implicit Solve Stage
    [UDGIMtmp,UHIMtmp] = hdg_solve(master,meshIM,app,UDGIMn,UHIMn,SDG);
    UDGIM = UDGIM + bAinv(istage)*UDGIMtmp;
    UHIM  = UHIM  + bAinv(istage)*UHIMtmp;
    
    % Transfer data to Explicit domain
    IMEXdata = getVhatIMEX(meshIM,UHIMtmp);
    
    % Add treatment of PDG for implicit domain (to be done)    
    if istage==1
        SPG = Ainv(1,1)*PDGIMn;
    elseif istage == 2
        SPG = (Ainv(2,1)+Ainv(2,2))*PDGIMn - Ainv(2,1)*PDGIMtmp;
    elseif istage == 3
        SPG = (Ainv(3,1)+Ainv(3,2)+Ainv(3,3))*PDGIMn - ...
              (Ainv(3,1)/Ainv(2,1))*((Ainv(2,1)+Ainv(2,2))*PDGIMn-SPG) - Ainv(3,2)*PDGIMtmp;
    end
    PDGIMtmp = (1/Ainv(istage,istage))*(UDGIMtmp(:,1,:)+SPG);
    PDGIM = PDGIM + bAinv(istage)*PDGIMtmp;
    
    
    % Compute new residual after Update
    ResEX(:,:,:,istage+1) = hdgres(master,meshEX,app,UDGEX);
    ResEXpdg(:,:,:,istage+1) = UDGEX(:,1,:);
    
end

% Update Explicit Solution at end of Time Step
ResSum = 0;
ResSumPdg = 0;
for j=1:nstage+1
    ResSum = ResSum + bhat(j)*ResEX(:,:,:,j);
    ResSumPdg = ResSumPdg + bhat(j)*ResEXpdg(:,:,:,j);
end
for i = 1:ne
    UDGEX(:,:,i) = UDGEXn(:,:,i) - dt*Minv(:,:,i)*ResSum(:,:,i); 
    PDGEX(:,:,i) = PDGEXn(:,:,i) + dt*ResSumPdg(:,:,i);
end


end



function [d,c,t] = dirkcoeff(q,p)
% q - number of stages
% p - order of accuracy

% fprintf('DIRK (%d,%d) \n', p,q);

if q == 1 && p == 1
    a = 1;
    b = 1;
    t = 1;
elseif q == 1 && p == 2
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
    a = 0.5;
    b = 1;
    t = 0.5;
elseif q == 2 && p == 2
    a = [1-0.5*sqrt(2), 0;
           0.5*sqrt(2), 1-0.5*sqrt(2)];
    b = [  0.5*sqrt(2), 1-0.5*sqrt(2)];
    t = [1-0.5*sqrt(2), 1];
elseif q == 2 && p == 3
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
    a = [0.5+0.5/sqrt(3), 0;
              -1/sqrt(3), 0.5 + 0.5/sqrt(3)];
    b = [            0.5, 0.5];
    t = [0.5+0.5/sqrt(3), 0.5-0.5/sqrt(3)];
elseif q == 3 && p == 3
    a1 = 0.4358665215;
    t1 = (1+a1)/2;
    b1 = -(6*a1^2-16*a1+1)/4;
    b2 = (6*a1^2-20*a1+5)/4;    
    a = [a1   ,  0, 0;
         t1-a1, a1, 0;
         b1   , b2, a1];
    b = [b1, b2, a1];
    t = [a1, t1, 1];    
elseif q == 3 && p == 4
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
    a1 = 2*cos(pi/18)/sqrt(3);
    a = [ 0.5*(1+a1),            0,          0;
             -0.5*a1,   0.5*(1+a1),          0;
                1+a1,    -(1+2*a1), 0.5*(1+a1)];
    b = [ 1/(6*a1^2), 1-1/(3*a1^2), 1/(6*a1^2)];
    t = [ 0.5*(1+a1),          0.5, 0.5*(1-a1)];
elseif q == 5 && p == 5
    a11 = (6-sqrt(6))/10;
    a21 = (-6+5*sqrt(6))/14;
    a31 = (888+607*sqrt(6))/2850; 
    a32 = (126-161*sqrt(6))/1425;
    a41 = (3153-3082*sqrt(6))/14250;
    a42 = (3213+1148*sqrt(6))/28500;
    a43 = (-267+88*sqrt(6))/500;
    a51 = (-32583+14638*sqrt(6))/71250;
    a52 = (-17199+364*sqrt(6))/142500;
    a53 = (1329-544*sqrt(6))/2500;
    a54 = (-96+131*sqrt(6))/625;
    b1  = 0;
    b2  = 0;
    b3  = (1/9);
    b4  = (16-sqrt(6))/36;
    b5  = (16+sqrt(6))/36;
    t1  = a11;
    t2  = (6+9*sqrt(6))/35;
    t3  = 1;
    t4  = (4-sqrt(6))/10;
    t5  = (4+sqrt(6))/10;
    a = [a11, 0, 0, 0, 0;
         a21, a11, 0, 0, 0;
         a31, a32, a11, 0, 0;
         a41, a42, a43, a11, 0;
         a51, a52, a53, a54, a11];
    b = [b1, b2, b3, b4, b5]; 
    t = [t1, t2, t3, t4, t5];            
else
    error('Invalid (q,p) combination');
end
d = inv(a);
c = b*d;

end


function [A,b,c] = ERKcoeff(q,p)
% function [D,b] = ERKcoeff(q,p)
% q - number of stages
% p - order of accuracy

if q == 1 && p == 1    
    A = 0;
    b = 1; 
    c = 1;
elseif q == 3 && p == 2
    alpha = 1 - (sqrt(2)/2);
    delta = -(2*sqrt(2))/3;
    
    A = [0,     0,       0;
         alpha, 0,       0;
         delta, 1-delta, 0];
    b = [0, 1-alpha, alpha]; 
    c = [0, alpha,    1];

elseif q == 3 && p == 3
    alpha = (3 + sqrt(3))/6;
    
    A = [0,       0,           0;
         alpha,   0,           0;
         alpha-1, 2*(1-alpha), 0];
    b = [0, 0.5,      0.5]; 
    c = [0, alpha,    1-alpha];    
    
elseif q == 4 && p == 3    
    A = [0,            0,            0,            0;
         0.4358665215, 0,            0,            0;
         0.3212788860, 0.3966543747, 0,            0;
         -0.105858296, 0.5529291479, 0.5529291479, 0];
     
     b = [0, 1.208496649, -0.644363171, 0.4358665215];
     c = [0, 0.4358665215, 0.7179332608, 1];
elseif q == 4 && p == 4    
    A = [0,            0,            0,            0;
         1/2, 0,            0,            0;
         0, 1/2, 0,            0;
         0, 0, 1, 0];
     
     b = [1/6, 1/3, 1/3, 1/6];
     c = [0, 1/2, 1/2, 1];     
else
    error('Invalid (q,p) combination');
end

%Dinv = A(2:end,1:end);
%Dinv = cat(1,Dinv,b);
%D = inv(Dinv);

end
