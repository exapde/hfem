function [err,erh] = faceerror2(UDG,UH,mesh,master)

%porder = mesh.porder;
perm = mesh.perm;
ne = mesh.ne;
nf = size(mesh.f,1);
nd = mesh.nd;
[npf,nfe] = size(perm);
ngf = master.ngf;
elcon = reshape(mesh.elcon,[npf nfe ne]);
shapft = squeeze(master.shapft(:,:,1));
%dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);

permgeom = master.permgeom; 
npfm = size(permgeom,1);
shapmf = squeeze(master.shapmf(:,:,1));
dshapmf  = reshape(permute(master.shapmf(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npfm]);

ncu = size(UH,1);
UH = reshape(UH,[ncu,npf,nf]);
err = zeros(ncu,1);
erh = err;
for i = 1:nf
    uh = reshape(UH(:,:,i),[ncu npf])';
    fi = mesh.f(i,end-1:end);
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
    i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
    j1 = elcon(:,i1,fi(1)) - (i-1)*npf;
%     pn = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));       
%     pg  = shapft*pn;        
%     dpg = dshapft*pn;    
    uhg = shapft*uh;
    
    pn = mesh.dgnodes(permgeom(:,i1),1:nd,fi(1));       
    pg  = shapmf*pn;        
    dpg = dshapmf*pn;           
    
    if nd==2
        jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
    elseif nd==3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
    end             
      
    %uhe = feval(func,pg);    
    if mesh.porder==0
        un = UDG(1,1:ncu,fi(1));     
    else
        un = UDG(perm(j1,i1),1:ncu,fi(1));     
    end
    uhe = shapft*un;
    an = abs(nlg(:,1)+nlg(:,2));
    aa = mean(nlg(:,1)+nlg(:,2));
    if aa<=0
        for j = 1:ncu        
            err(j) = err(j) + (master.gwfc.*jac)'*(an.*(uhg(:,j)-uhe(:,j)).^2);  
        end    
    end
    
    if  (mesh.f(i,end)==-2 || mesh.f(i,end)==-3)        
        if mean(pg(:,2))<=0.8
            uhe = 0;
        else
            uhe = 1;
        end
        for j = 1:ncu
            erh(j) = erh(j) + (master.gwfc.*jac)'*(an.*(uhg(:,j)-uhe).^2);  
        end            
    end
end
err  = sqrt(err);
erh  = sqrt(erh);


