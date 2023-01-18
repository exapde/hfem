function [UDG,UH,Nit]=hdg_solve(master,mesh,app,UDG,UH,SH,NitMax)
%HDG_SOLVE Solve using the HDG method and Newton iteration
%   [UH,QH,PH,UHAT] = HDG_SOLVE(MASTER,MESH,UDG,UH,SH,APP)
%
%      MASTER:                  Master structure
%      MESH:                    Mesh structure
%      UH(NPL,NC,NE):           Vector of unknowns (initial guess)
%      QH(NPL,NC,2,NE):         Vector of gradients of U (initial guess)
%      PH(NPL,1,NE):            Vector of pressure (initial guess)
%      UHAT(NC,3*NPS,NF):       Vector of U_hat's (initial guess)
%      APP:                     Application structure
%      SH:                      Source term independent of U (for time dependent mode only)
%
%      UH(NPL,NC,NE):           Vector of unknowns
%      QH(NPL,NC,2,NE):         Vector of gradients of U
%      PH(NPL,1,NE):            Vector of pressure
%      UHAT(NCH,3*NPS,NF):      Vector of U_hat's
%
%      NPL:                     Number of DG nodes within an element
%      NC:                      Number of conservation equations solved (components)
%      NPS:                     Number of HDG nodes per edge (porder+1)
%      NE:                      Number of elements
%      NF:                      Number of faces
%

if nargin < 7
    NitMax=8;
end

npf = master.npf;
npv = master.npv;
ne  = mesh.ne;
%nf  = mesh.nf;
nfe = size(master.perm,2);
nc  = app.nc;
ncu = app.ncu;
%nd  = app.nd;
nsiz = mesh.nsiz;

if isfield(app,'fbou') == 0
    app.fbou = 'fbou';
end
if isfield(app,'source') == 0
    app.source = 'source';
end
if isfield(app,'denseblock') == 0
    app.denseblock = 0;
end
if isfield(app,'linear') == 0
    app.linear = 0;
end
if isfield(app,'adjoint') == 0
    app.adjoint = 0;
end
if isfield(app,'getdqdg') == 0
    app.getdqdg = 1;
end

if app.getdqdg == 0
    nco = ncu;
else
    nco = nc;
end

% Imposing the source Term SH
if isempty(SH)
    SH = zeros(npv,nc,ne);
end
if app.wave==0
    SH(:,ncu+1:end,:)=0;
else
    if app.flg_p
        SH(:,ncu+1,:)=0;
    end
end

if app.denseblock ~= 0
    f2f = mkf2f(mesh.f, mesh.t2f);
end

if min(master.perm(:))==0
    master.perm = master.perm + 1;
end
if min(master.permgeom(:))==0
    master.permgeom = master.permgeom + 1;
end

tic
[K,F,DUDG,DUDG_DUH] = hdg_assemble(master,mesh,app,UDG,UH,SH);  
toc

% For Adjoint problems only
if app.adjoint == 1 
    if app.denseblock == 0
        UHE  = full(reshape(K\F,ncu,nsiz));
    else
        UHE = gmres(K, F, f2f);
        UHE = full(reshape(UHE,ncu,nsiz));    
    end
    if isempty(DUDG) == 0
        UHE = reshape(full(UHE(:,mesh.elcon(:))),ncu,nfe*npf,ne);    
        UDG = DUDG + reshape(mapContractK(DUDG_DUH,UHE,[1 2],[3 4],5,[1 2],[],3),[npv nco ne]);
    else        
        UDG = hdg_local(master,mesh,app,UDG,UH,SH,UHE);        
    end
    UH = UHE;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      NEWTON-RAPHSON LOOP                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
it   = 0;
duh  = 1e6;
NewtonTol = 1e-8;

while duh > NewtonTol && it < NitMax
    
    % Solving for traces increments DUH via the Schurr complement
    if app.denseblock==0
        %DUH  = sparsesolve(K, F, mesh.f, mesh.f2f, 0, app.hybrid);
        %DUH = full(reshape(DUH,ncu,nsiz));
        DUH = full(reshape(K\F,ncu,nsiz));
    else
        DUH = gmres(K, F, f2f);
        DUH = full(reshape(DUH,ncu,nsiz));    
    end    
    
    % Computing the variables increments DUDGT
    if isempty(DUDG)==0
        % Reorganize DUH vector using connectivities
        DUHE = reshape(full(DUH(:,mesh.elcon(:))),ncu,nfe*npf,ne);
        % Update increment DUDG from traces increments
        DUDGT = DUDG + reshape(mapContractK(DUDG_DUH,DUHE,[1 2],[3 4],5,[1 2],[],3),[npv nco ne]);
    
    else % compute DUDGT by solving the local problems
        DUDGT = hdg_local(master,mesh,app,UDG,UH,SH,DUH);
    end
    
    % Linear Flag, currently unused for linear elasticity
    if app.linear == 1
        UDG(:,1:nco,:)  = UDG(:,1:nco,:) + DUDGT;                     
        UH   = UH + DUH;   
        if app.getdqdg==0
            QDG = getq(master, mesh, UDG, UH, SH, app.fc_q);
            UDG(:,ncu+1:end,:)=QDG;
        end        
        return; 
    end        
    
    it = it + 1;    
    %fprintf('Newton iteration :  %d\n', it);
    
    % Find stepsize for damped Newton    
    UDG0 = UDG;
    UH0  = UH;    
    duh0 = norm(F(:));   
    alfa = 1;
    while 1
        UDG(:,1:nco,:)  = UDG0(:,1:nco,:) + alfa*DUDGT;                     
        UH   = UH0 + alfa*DUH;   
        
        if app.getdqdg==0
            QDG = getq(master, mesh, UDG, UH, SH, app.fc_q);
            UDG(:,ncu+1:end,:)=QDG;
        end        
        
        [K,F,DUDG,DUDG_DUH] = hdg_assemble(master,mesh,app,UDG,UH,SH);        
        duh  = norm(F(:));                 
        
        %if it==1
        %    break;
        %end
        
        if ((duh>duh0) || isnan(duh))
            alfa=alfa/2;
            %fprintf(' alpha = %f\n',alfa);
            if alfa<5e-1, break; end
            %if alfa<1e-2, break; end
        else
            break;
        end
    end
       
    if duh>1.e3
        it = NitMax;
    end
    if nargout > 2
        Nit = it;
    end
    fprintf('Iteration: %d,   Old residual: %e,   New residual: %e    %e\n', [it duh0 duh alfa]);
    %fprintf('Iteration: %d,   New residual: %e,    alpha: %e\n', [it duh alfa]);
end
fprintf('Newton iterations :  %d\n', it); % a supprimer !!!
