function [y,u,v,w] = uhplot(mesh,master,UH,UDG)

npf = size(mesh.permgeom,1);
nfe = size(mesh.permgeom,2);
ne  = size(mesh.dgnodes,3);
nd  = mesh.nd;
nf = mesh.nf;

%UH  = reshape(UH(mesh.elcon),[nqf nfe ne]);
UH  = reshape(UH,[],nf);
%elcon = reshape(mesh.elcon,[npf nfe ne]);

xi     = linspace(0,1,40)';
[~,~,plocfc] = masternodes(max(1,master.porder),mesh.nd,mesh.elemtype,mesh.nodetype);   
shapmf = mkshape(max(1,master.porder),plocfc,xi,1);
shapmf = shapmf(:,:,1)';

x = [];
u = [];
v = [];
w = [];
for k = 1:nf
    uh = UH(:,k);
    if mesh.porder==0
        uh = [uh; uh];
    end
    uh = shapmf*uh;
    
    fi = mesh.f(k,end-1); % obtain two elements sharing the same face i  
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces     
    i1 = find(kf(1,:)==k);  % obtain the index of face i in the 1st element                     
    dg = mesh.dgnodes(mesh.permgeom(:,i1),1:nd,fi(1));                                                                 
    dg = shapmf*reshape(dg,[npf nd]);
    udg = UDG(mesh.perm(:,i1),1,fi(1));   
    if mesh.porder==0
        udg = [udg; udg];
    end
    udg = shapmf*reshape(udg,[npf 1]);    
    u1 = uh + 0.5*(udg-uh);
    u2 = uh - 0.5*(udg-uh);
    %dg = reshape(dg,[],nd);
    if max(abs(dg(:,1)+dg(:,2)-1))<1e-3
        %plot3(dg(:,1),dg(:,2),uh(:),'-k');
        x = [x; dg(:,1) dg(:,2)];
        u = [u; u1(:)];
        v = [v; u2(:)];
        w = [w; uh(:)];
    end
end
y = sqrt((x(:,1)-0).^2 + (x(:,2)-1).^2);





