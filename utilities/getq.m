function QDG = getq(mesh,mastersubgrid,UDG,UHAT)

if iscell(mastersubgrid)==0
    QDG = calq(mastersubgrid,mesh.dgnodes,UHAT,UDG);
    return;
end

for j = 1:length(mastersubgrid) % for each subgrid
    % obtain element indices for subgrid j
    ind = find(mesh.subgrids==j);    
    if isempty(ind)==0        

        % obtain geometry nodes
        dgnodes = mesh.dgnodes(:,:,ind);                
        npm = size(dgnodes,1);
        nd  = size(dgnodes,2);
        nj  = size(dgnodes,3);
        
        % obtain subgrid geometry nodes
        geomnodes = mastersubgrid{j}.geomnodes;        
        ng = size(geomnodes,3);

        % compute shape functions at the subgrid DG nodes
        shp = zeros(npm,npm,ng);                
        for k=1:ng
            tmp = mkshape(mesh.porder,mesh.plocal,geomnodes(:,:,k),mesh.elemtype);
            shp(:,:,k) = tmp(:,:,1)';
        end      

        % compute subgrid DG nodes in the physical space
        dgx = zeros(npm,nd,nj*ng);
        for i=1:nj % for each superelement                           
            for k=1:ng % for each subgrid element                    
                dgx(:,:,(i-1)*ng+k) = shp(:,:,k)*dgnodes(:,:,i);                   
            end
        end        
        
        ncu = size(UHAT{j},1);
        npf  = size(mastersubgrid{j}.perm,1);
        nfe  = size(mastersubgrid{j}.perm,2);
        uhat = reshape(UHAT{j}(:,mastersubgrid{j}.elcon,:),[ncu npf*nfe nj*ng]);    
        QDG{j} = calq(mastersubgrid{j},dgx,uhat,UDG{j});
    end
end
        
function Q = calq(master,dgnodes,UHAT,UDG)

nd   = master.nd;
npv  = size(UDG,1);
ne   = size(UDG,3);
npf  = size(master.perm,1);
nfe  = size(master.perm,2);
ncu  = numel(UHAT)/(npf*nfe*ne);

[~, minusJeInvDetJe, detJe] = volgeom(master,dgnodes);
[~, nlgf, detJf] = facegeom(master,dgnodes);
[MiC, MiE] = qint(master,detJe,minusJeInvDetJe,detJf,nlgf);
MiC = permute(MiC,[1 3 2 4]);
MiE = permute(MiE,[1 3 2 4]);

Q = zeros(npv,ncu*nd,ne);
for k=1:ne                          
    u    = UDG(:,1:ncu,k);
    uhat = reshape(UHAT(:,:,k),[ncu npf*nfe]);             
    MiCt = reshape(MiC(:,:,:,k), [npv*nd npv]);
    MiEt = reshape(MiE(:,:,:,k), [npv*nd npf*nfe]);
    MiCu = permute(reshape(MiCt*u,[npv nd ncu]), [1 3 2]);
    MiEu = permute(reshape(MiEt*(uhat'),[npv nd ncu]), [1 3 2]);
    Q(:,:,k) = reshape(MiEu,[npv ncu*nd]) - reshape(MiCu,[npv ncu*nd]);
end

