function plotudg(mesh,master,UDG)

npf = size(mesh.permgeom,1);
nd  = mesh.nd;
nf = mesh.nf;

xi     = linspace(0,1,40)';
[~,~,plocfc] = masternodes(max(1,master.porder),mesh.nd,mesh.elemtype,mesh.nodetype);   
shapmf = mkshape(max(1,master.porder),plocfc,xi,1);
shapmf = shapmf(:,:,1)';

figure(1); clf;
for k = 1:nf    
    fi = mesh.f(k,end-1:end); % obtain two elements sharing the same face i  
    kf = mesh.t2f(fi,:);    % obtain neighboring faces     
    i1 = find(kf(1,:)==k);  % obtain the index of face i in the 1st element   
    i2 = kf(2,:)==k;  % obtain the index of face i in the 1st element                     
    dg = mesh.dgnodes(mesh.permgeom(:,i1),1:nd,fi(1));                                                                 
    dg = shapmf*reshape(dg,[npf nd]);
    udg1 = UDG(mesh.perm(:,i1),1,fi(1));   
    udg2 = UDG(mesh.perm(:,i2),1,fi(2));   
    udg1 = shapmf*reshape(udg1,[npf 1]);    
    udg2 = shapmf*reshape(udg2,[npf 1]);    
    udg = 0.5*(udg1+udg2);    
    if max(abs(dg(:,1)))<1e-10
        plot(dg,udg);
    end
end






