function [mastersubgrid,UDG,UHAT] = initm(mesh,ui,ncu)

nc  = length(ui);
for i = 1:size(mesh.gridtype,1)
    mastersubgrid{i} = mkmastersubgrid(mesh.porder,mesh.gridtype(i,1),max(2*mesh.gridtype(i,1),1),mesh.gridtype(i,2),mesh.nd,mesh.elemtype,mesh.nodetype);        
    npv  = size(mastersubgrid{i}.plocvl,1);
    npf  = size(mastersubgrid{i}.perm,1);        
    nf   = max(mastersubgrid{i}.elcon(:))/(npf);
    nes  = size(mastersubgrid{i}.t,1);    
    ind  = find(mesh.subgrids==i);
    ne   = length(ind);
    UDG{i} = zeros(npv,nc,nes*ne);
    UHAT{i} = zeros(ncu,npf*nf,ne);    
    for j=1:nc
        UDG{i}(:,j,:) = ui(j);
    end    
    for k = 1:ne
        UHAT{i}(:,:,k) = inituhat(mastersubgrid{i},mastersubgrid{i}.elcon,UDG{i}(:,:,(k-1)*nes+1:k*nes),ncu);  
    end    
end

