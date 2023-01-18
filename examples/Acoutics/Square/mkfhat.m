function [UHAT,QHAT] = mkfhat(app,mesh,UDG,tau)

time = app.time;
kx  = app.k(1);
ky  = app.k(2);
kk  = sqrt(kx^2+ky^2);

shapfc = mkshape(mesh.porder,mesh.plocfc,mesh.plocfc,mesh.elemtype);
dshapft  = shapfc(:,:,2)';

nd = mesh.nd;
[npf,nfe] = size(mesh.perm);
[~,~,ne] = size(UDG);
perm = mesh.perm;
t2t = mesh.t2t; %mkt2t(mesh.t,mesh.elemtype);
UHAT = zeros(npf,nfe,ne);
QHAT = zeros(npf,nfe,ne);
for i = 1:ne    
    e1 = i;    
    for j = 1:nfe % for each face of element i        
        e2 = t2t(i,j);
        u1 = UDG(perm(:,j),:,e1);        
        pb = mesh.dgnodes(perm(:,j),:,e1);
        dpg = reshape(dshapft*pb,[npf nd]);     
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nl   = [dpg(:,2)./jac,-dpg(:,1)./jac];
        if e2>0                                                
            k  = mesh.t2f(e2,:)==mesh.t2f(e1,j);  % obtain the index of face i in the second element        
            u2 = UDG(perm(end:-1:1,k),:,e2);         
            qn = sum((u1(:,2:end)-u2(:,2:end)).*nl,2);
            UHAT(:,j,i) = 0.5*(u1(:,1)+u2(:,1)) - (0.5/tau)*qn;            
            qn = sum(u1(:,2:end).*nl,2);
            QHAT(:,j,i) = qn - tau*(u1(:,1) - UHAT(:,j,i));
        else % boundary conditions            
            UHAT(:,j,i) = 0;
            qn = sum(u1(:,2:end).*nl,2);
            QHAT(:,j,i) = qn - tau*(u1(:,1) - UHAT(:,j,i));
            
%             % absorbing condition
%             if mesh.f(abs(mesh.t2f(i,j)),end)==-1 || mesh.f(abs(mesh.t2f(i,j)),end)==-2 || mesh.f(abs(mesh.t2f(i,j)),end)==-3 || mesh.f(abs(mesh.t2f(i,j)),end)==-4
%                 qn = sum(u1(:,2:end).*nl,2);
%                 UHAT(:,j,i) = (1/(1+tau))*(tau*u1(:,1) - qn);                
%                 QHAT(:,j,i) = qn - tau*(u1(:,1) - UHAT(:,j,i));                
%             end
%             
%             % reflection condition
%             if mesh.f(abs(mesh.t2f(i,j)),end)==-5
%                 qn = sum(u1(:,2:end).*nl,2);                
%                 q1 = (kx)*cos(kx*pb(:,1)+ky*pb(:,2) - kk*time); 
%                 q2 = (ky)*cos(kx*pb(:,1)+ky*pb(:,2) - kk*time); 
%                 gn = q1.*nl(:,1) + q2.*nl(:,2);         
%                 UHAT(:,j,i) = (1/(tau))*(tau*u1(:,1) - qn - gn);         
%                 QHAT(:,j,i) = qn - tau*(u1(:,1) - UHAT(:,j,i));                
%             end
        end                 
    end
end
UHAT = reshape(UHAT,[npf*nfe,ne]);
QHAT = reshape(QHAT,[npf*nfe,ne]);

