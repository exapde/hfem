
function [Fx,Fy] = getliftdrag(master,mesh,app,UDG,UH,wid,ix,iy)

nd  = master.nd;
ngf = master.ngf;
npf = master.npf;
ne  = size(mesh.t,1);
nfe = size(master.perm,2);

elcon  = mesh.elcon;
perm   = master.perm(:,:,1);
shapft = squeeze(master.shapft(:,:,1));
dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
gwfc   = master.gwfc;

Fx=0;
Fy=0;

for i=1:ne
    dg = mesh.dgnodes(:,1:nd,i);       
    uh = UH(:,elcon(:,i))';          % uh(3*nps,nc)
    udg = UDG(:,:,i);
    
    fc=mesh.f(abs(mesh.t2f(i,:)),end);  % fc=mesh.f(abs(mesh.t2f(i,:)),4);
    bf  = (fc<0);
    
    for is = 1:nfe
        if bf(is)==1 % boundary face
%             ib = app.bcm(-fc(is));
            if fc(is) == -wid %ib == wid
                pn    = dg(perm(:,is),:);
                pg    = shapft*pn;
                udgg  = shapft*udg(perm(:,is),:);
                uhg   = shapft*uh((is-1)*npf+1:is*npf,:);                                    
                dpg   = dshapft*pn;
                dpg   = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 2]);                
                
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
                
                [fh] = fhat(nlg,pg,udgg,uhg,app.arg,0);
                Fx = Fx + sum(fh(:,ix).*jac(:).*gwfc(:));                 
                Fy = Fy + sum(fh(:,iy).*jac(:).*gwfc(:));                 
            end
        end   
    end
end
