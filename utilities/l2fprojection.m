function UH = l2fprojection(mesh,master,func,param,time,nc)

nd  = master.nd;
nf = mesh.nf;
ne = mesh.ne;
npf = master.npf;
ngf = master.ngf;
nfe  = size(master.perm,2);

elcon = reshape(mesh.elcon,[npf nfe ne]);
perm   = master.perm(:,:,1);
shapft = squeeze(master.shapft(:,:,1));
dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);

UH = zeros(npf,nc,nf);
for i = 1:nf
    fi = mesh.f(i,end-1:end);    
    kf = mesh.t2f(fi(1),:);  
    i1 = find(kf(1,:)==i);  
    j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
    pn = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));        
    
    pg    = shapft*pn;
    dpg   = dshapft*pn;    
    if nd==2
        jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);        
    elseif nd==3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);        
    end             
                                
    fg = feval(func,pg,param,time);        
    M = reshape(master.shapfgdotshapfc(:,:,1)*jac,[npf npf]);        
    for j = 1:nc
        F = master.shapfg(:,:,1)*(jac.*fg(:,j));
        UH(:,j,i) = M\F;        
    end    
%     [UH(:,:,i) feval(func,pn,param,time)]
%     pause
end

UH = permute(UH, [2 1 3]);
UH = reshape(UH, [nc npf*nf]);


