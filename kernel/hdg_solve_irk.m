function [UDGn,UHn,U] = hdg_solve_irk(master,mesh,app,UDG,UH,U0,time,dt,nstage,torder)
%HDG_SOLVE_DIRK Solve using the HDG method and Newton's iteraion - DIRK timestepping
%   [UH,QH,UHAT] = HDG_SOLVE_DIRK(MASTER,MESH,UH,QH,UHAT,APP,TIME,DT,NSTAGE,TORDER)
%
%      MASTER:                  Master structure
%      MESH:                    Mesh structure
%      UH(NPL,NC,NE):           Vector of unknowns (initial guess)
%      QH(NPL,NC,2,NE):         Vector of gradients of U (initial guess)
%      UHAT(NC,3*NPS,NF):       Vector of U_hat's (initial guess)
%      APP:                     Application structure
%      TIME:                    Time
%      DT:                      Timestep
%      NSTAGE:                  Number of stages
%      TORDER:                  Order of accuracy
%
%      UH(NPL,NC,NE):           Vector of unknowns
%      QH(NPL,NC,2,NE):         Vector of gradients of U
%      UHAT(NC,3*NPS,NF):       Vector of U_hat's
%
%      NPL:                     Number of DG nodes within an element
%      NC:                      Number of conservation equations solved (components)
%      NPS:                     Number of HDG nodes per edge (porder+1)
%      NE:                      Number of elements
%      NF:                      Number of faces
%
%  Ref: Diagonally Implicit Runge Kutta Methods for Stiff ODE's, 
%  by Roger Alexander, SINUM, vol. 14, no. 6, 1997.
%

ncu = app.ncu;
nch0 = ncu/nstage;
nc0 = app.nc/nstage;
[D,t,c] = irkcoeff(nstage,torder);
app.arg{2}.D = (1/dt)*D;
app.arg{2}.t = t;
app.arg{2}.dt = dt;
app.arg{2}.Ns = nstage;
app.arg{2}.nch0 = ncu/nstage;
app.arg{2}.nd = app.nd;
app.fc_u = 0;
app.fc_q = 1;
app.time = time;

% initial solution
% nch0 = ncu/nstage;
% i = nstage;
% uvars0 = (i-1)*nch0+1:i*nch0; 
% U0 = UDG(:,uvars0,:);

% compute source term due to U0
SDG = 0*UDG;
for i=1:nstage
    % subset of U variables describing this snapshot
    uvars = (i-1)*nch0+1:i*nch0; 
    
    for j=1:nch0
        SDG(:,uvars(j),:) = sum(D(i,:)/dt)*U0(:,j,:);                
    end        
end

% for i=1:Ns
%     % subset of U variables describing this snapshot
%     uvars = (i-1)*nch0+1:i*nch0; 
%     
%     % subset of UDG variables describing this snapshot
%     udgvars = zeros(nch0,nd1);
%     for j=1:nch0; udgvars(j,:) = (i-1)*nch0+j:nch:nc; end
%     udgvars= sort(udgvars(:))';
%     
%     % time coordinate of this snapshot    
%     timei = time + t(i)*dt;    
%     
%     % source term S
%     [sr(:,uvars),sr_udg(:,uvars,udgvars)] = source(p,udg(:,udgvars),param{1},timei);
%     
%     % source term -D*U 
%     for j=1:nch0
%         sr(:,uvars(j)) = sr(:,uvars(j)) - udg(:,j:nch0:nch)*(D(i,:)');
%         sr_udg(:,uvars(j),j:nch0:nch) = sr_udg(:,uvars(j),j:nch0:nch) ...
%           - permute(repmat(D(i,:),ng,1).*ones(size(udg(:,j:nch0:nch))),[1 3 2]);           
%     end        
% end

% for i=1:Ns
%     % subset of U variables describing this snapshot
%     uvars = (i-1)*nch0+1:i*nch0; 
%     
%     % subset of UDG variables describing this snapshot
%     if nd1==1
%         udgvars = uvars;
%     else
%         udgvars = zeros(nch0,nd1);
%         for j=1:nch0; udgvars(j,:) = (i-1)*nch0+j:nch:nc; end
%         udgvars= sort(udgvars(:))';
%     end
%     
%     % time coordinate of this snapshot    
%     timei = time + t(i)*dt;
%     
%     % Assemble fluxes
%     [f(:,uvars,:),f_udg(:,uvars,:,udgvars)] = flux(p,udg(:,udgvars),param{1},timei);
% end

% HDG solver
[UDGn,UHn] = hdg_solve(master,mesh,app,UDG,UH,SDG);

nc = app.nc;
nch = app.nch;
nd1 = nc/nch;
U = (1-sum(c))*U0;
for i=1:nstage
    
    udgvars = zeros(nch0,nd1);
    for j=1:nch0; udgvars(j,:) = (i-1)*nch0+j:nch:nc; end
    udgvars= sort(udgvars(:))';
    
    U = U + c(i)*UDGn(:,udgvars,:);
end

function [d,t,c] = irkcoeff(q,p)
% q - number of stages
% p - order of accuracy

if q == 1 && p == 1 % BDF
    a = 1;
    b = 1;
    t = 1;
elseif q == 1 && p == 2 % Gauss
    a = 0.5;
    b = 1;
    t = 0.5;    
elseif q == 2 && p == 2 % Lobatto IIIC
    a = [1/2, -1/2;
         1/2,  1/2];
    b = [1/2 1/2];
    t = [0, 1];            
elseif q == 2 && p == 3 % RADAU IIA    
    a = [5/12, -1/12;
         3/4,   1/4];
    b = [3/4 1/4];
    t = [1/3, 1];        
elseif q == 2 && p == 4 % Gauss 
    a = [1/4, 1/4-sqrt(3)/6;
         1/4+sqrt(3)/6,   1/4];
    b = [1/2 1/2];
    t = [1/2-sqrt(3)/6, 1/2+sqrt(3)/6];                
elseif q == 3 && p == 4  % Lobatto IIIC        
    a = [1/6, -1/3,  1/6;
         1/6, 5/12, -1/12;
         1/6,  2/3,  1/6];
    b = [1/6 2/3 1/6];
    t = [0 1/2 1];        
elseif q == 3 && p == 5  % RADAU IIA    
    a = [(88-7*sqrt(6))/360  (296-169*sqrt(6))/1800 -(2-3*sqrt(6))/225;
         (296+169*sqrt(6))/1800  (88+7*sqrt(6))/360 -(2+3*sqrt(6))/225;
         (16-sqrt(6))/36 (16+sqrt(6))/36 1/9];
    b = [(16-sqrt(6))/36 (16+sqrt(6))/36 1/9];
    t = [(4-sqrt(6))/10 (4+sqrt(6))/10 1];       
elseif q == 3 && p == 6  % Gauss
    a = [5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
         5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24;
         5/36+sqrt(15)/30,  2/9+sqrt(15)/15,  5/36];
    b = [5/18 8/18 5/18];
    t = [1/2-sqrt(15)/10, 1/2, 1/2+sqrt(15)/10];    
else
    error('Invalid (q,p) combination');
end
d = inv(a);
c = b*d;

