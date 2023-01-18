function [AE,FE,DUDG,DUDG_DUH] = hdg_assemble(master,mesh,app,UDG,UH,SH)


ne = mesh.ne;
npv = master.npv;
npf = master.npf;
nfe = size(master.perm,2);
ncu = app.ncu;
nc  = app.nc;
nb  = npf*ncu;

% mesh partioning by metis
%epart = meshpart(mesh.t,np);
%np = max(1,matlabpool('size'));
np = 1;
nm = floor(ne/np);
ns = nm*ones(np,1);
for j=1:ne-nm*np
    ns(j) = ns(j)+1;
end
epart = zeros(ne,1);
for i=2:np
    ind = sum(ns(1:i-1))+1:sum(ns(1:i));
    epart(ind) = i-1;
end
% make slice variables for parfor
for i = 1:np    
    ind = find(epart==i-1);  
    masterp{i} = master;
    appp{i} = app;                  
    udg{i} = UDG(:,:,ind);
    sh{i}  = SH(:,:,ind);
    uh{i}  = reshape(UH(:,mesh.elcon(:,ind)),[ncu npf*nfe length(ind)]);
    dgnodes{i} = mesh.dgnodes(:,:,ind);
    bf{i} = mesh.bf(:,ind);
end
clear UDG SH UH 

% compute elemental matrices in parallel
for i = 1:np 
    [ae{i},fe{i},dudg{i},dudg_duh{i}] = hdg_elemental(masterp{i},appp{i},dgnodes{i},bf{i},udg{i},uh{i},sh{i});
end


if app.getdqdg == 1
    nco = nc;
elseif app.getdqdg == 0
    nco = ncu;
else
    nco = 0;
end
% assemble into single matrices
AE  = zeros(ncu,npf*nfe,ncu,npf*nfe,ne);
FE  = zeros(ncu,npf*nfe,ne);
DUDG_DUH = zeros(npv,nco,ncu,npf*nfe,ne);
DUDG = zeros(npv,nco,ne);
for i=1:np
    ind = find(epart==i-1);    
    
    DUDG(:,:,ind) = dudg{i};
    DUDG_DUH(:,:,:,:,ind) = dudg_duh{i};
    AE(:,:,:,:,ind) = ae{i};
    FE(:,:,ind) = fe{i};
    
    dudg{i} = [];
    dudg_duh{i} = [];
    ae{i} = [];
    fe{i} = [];
end

if app.denseblock == 1
    % assemble the matrix system in dense format
    [AE,FE] = mkDenseBlockSystem(AE, FE, mesh.f, mesh.t2f, mesh.elcon);    
else
    % assemble the matrix system in sparse format         
    il = zeros(nfe*npf*ncu,nfe*npf*ncu,ne);
    jl = zeros(nfe*npf*ncu,nfe*npf*ncu,ne);
    for i=1:ne
        con = repmat((mesh.elcon(:,i)'-1)*ncu,ncu,1)+repmat((1:ncu)',1,nfe*npf);
        con = reshape(con,nfe*npf*ncu,1);    
        il(:,:,i) = repmat(con ,1,nfe*nb);
        jl(:,:,i) = repmat(con',nfe*nb,1);        
    end

    AE = sparse(reshape(il,nfe*nfe*nb*nb*ne,1),reshape(jl,nfe*nfe*nb*nb*ne,1),reshape(AE,nfe*nfe*nb*nb*ne,1));        
    FE = sparse(reshape(il(:,1,:),nfe*nb*ne,1),ones(nfe*nb*ne,1),reshape(FE,nfe*nb*ne,1));          
end