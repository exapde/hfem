function DUDGT = hdg_local(master,mesh,app,UDG,UH,SH,DUH)

ne = mesh.ne;
npv = master.npv;
npf = master.npf;
nfe = size(master.perm,2);
nc  = app.nc;
ncu = app.ncu;
np  = app.np;

% mesh partioning by metis
%epart = meshpart(mesh.t,np);
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
    duh{i}  = reshape(DUH(:,mesh.elcon(:,ind)),[ncu npf*nfe length(ind)]);
    dgnodes{i} = mesh.dgnodes(:,:,ind);
    bf{i} = mesh.bf(:,ind);
end
clear UDG SH UH DUH

parfor i = 1:np 
    dudgt{i} = localsolve(masterp{i},appp{i},dgnodes{i},bf{i},udg{i},uh{i},sh{i},duh{i});
end

DUDGT = zeros(npv,nc,ne);
for i=1:np
    ind = (epart==i-1);        
    DUDGT(:,:,ind) = dudgt{i};    
    dudgt{i} = [];    
end


function DUDGT = localsolve(master,app,dgnodes,bf,UDG,UH,SH,DUH)

npv = size(UDG,1);                
ne  = size(dgnodes,3);
nc  = app.nc;
ns  = 200;

nb = ceil(ne/ns);          
nk = 1:ns:ne;
nm = [nk(1:end); [nk(2:end)-1,ne]];    

UDG = permute(UDG,[1 3 2]);
SH  = permute(SH,[1 3 2]);
dgnodes = permute(dgnodes,[1 3 2]);
UH = permute(UH,[2 3 1]);

DUDGT = zeros(npv,nc,ne);
for j=1:nb
    id = nm(1,j):nm(2,j);
    
    % compute volume integrals
    [Ru, Rq, BD, M, C, L, Q, Ju, Jq] = volint(master,app,dgnodes(:,id,:),UDG(:,id,:),SH(:,id,:));
    
    % compute face integrals
    [Ru, Rq, Rh, BD, F, E, GK, H, Ju, Jq, Jh] = faceint(master,app,bf(:,id),dgnodes(:,id,:),UDG(:,id,:),UH(:,id,:),Ru,Rq,BD,Ju,Jq);
    
    % perform schur compliment to obtain the elemental matrices and vectors
    if app.adjoint == 0
        [dudg, dudg_duh] = schur_primal(M, C, E, L, Q, BD, F, GK, H, Rq, Ru, Rh);
    else
        [dudg, dudg_duh] = schur_adjoint(M, C, E, L, Q, BD, F, GK, H, Jq, Ju, Jh);
    end
        
    DUDGT(:,:,id) = dudg + reshape(mapContractK(dudg_duh,DUH(:,:,id),[1 2],[3 4],5,[1 2],[],3),[npv nc length(id)]);
end    
