function [M,Minv]=massinv(master, mesh)

nd  = size(mesh.p,2);
ne  = size(mesh.t,1);
npv = size(mesh.dgnodes,1);
ngv = master.ngv;

dshapvt(:,:,1) = master.shapvl(:,:,2)';
dshapvt(:,:,2) = master.shapvl(:,:,3)';
dshapvt = reshape(permute(dshapvt,[1 3 2]),[ngv*nd npv]);

shapvgdotshapvl  = zeros(npv*npv,ngv,nd+1);      
for d=1:nd+1    
    shapvg = master.shapvl(:,:,d)*diag(master.gwvl);    
    for ii=1:npv
        for jj = 1:npv
            shapvgdotshapvl((ii-1)*npv+jj,:,d) = shapvg(jj,:).*master.shapvl(ii,:,1);                    
        end
    end            
end

M = zeros(npv,npv,ne);
Minv = zeros(npv,npv,ne);
nblk(1,:) = 1:100:ne;
nblk(2,:) = [nblk(1,2:end) ne];
for n = 1:size(nblk,2)
    e1 = nblk(1,n);
    e2 = nblk(2,n);                         
    ns = (e2-e1+1);
    
    pn = reshape(mesh.dgnodes(:,1:nd,e1:e2),[npv nd*ns]);    
    Jg = reshape(dshapvt*pn,[ngv nd nd ns]);  
    Jg = permute(Jg,[1 4 2 3]);
    jac = Jg(:,:,1,1).*Jg(:,:,2,2) - Jg(:,:,1,2).*Jg(:,:,2,1);    
        
    tm = shapvgdotshapvl(:,:,1)*reshape(jac,[ngv ns]);
    M(:,:,e1:e2) = reshape(tm,[npv npv ns]);
    for i = e1:1:e2
        Minv(:,:,i) = inv(M(:,:,i));
    end
end
