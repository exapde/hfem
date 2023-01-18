function [dgnodes,u,uint] = surfacedata(mesh,master,UDG,ib)

perm = mesh.perm;
ne = mesh.ne;
nd = mesh.nd;
[npf,nfe] = size(perm);
elcon = reshape(mesh.elcon,[npf nfe ne]);

in = find(mesh.f(:,end)==ib);
if isempty(in)
    error('Boundary is invalid.');
end
ns = length(in); 

shapft = squeeze(master.shapft(:,:,1));
ngf = size(shapft,1);
dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);

uint = 0;
dgnodes = zeros(ngf,nd,ns);
u = zeros(ngf,ns);
for j = 1:ns
    i = in(j);
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
    i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
    j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
    p = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));            
    dgnodes(:,:,j) = shapft*p;
    u(:,j) = shapft*UDG(perm(j1,i1),1,fi(1));       
    
    dpg   = dshapft*p;
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
    
    uint = uint + sum(abs(jac).*master.gwfc.*u(:,j));
end


