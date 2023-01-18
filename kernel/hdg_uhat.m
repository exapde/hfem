function UH = hdg_uhat(mesh,master,app,UDG)

nd = mesh.nd;
%nc = app.nc;
ncu = app.ncu;   % number of components of UH
nf  = size(mesh.f,1);    % number of faces
npf = size(mesh.perm,1); % number of points per face
ngf = master.ngf;

%ubou = str2func(app.ubou);
% shapfc = mkshape(mesh.porder,master.plocfc,master.plocfc,mesh.elemtype);
% dshapft  = reshape(permute(shapfc(:,:,2:nd),[1 3 2]),[npf*(nd-1) npf]);

% Shap functions 
perm            = master.perm(:,:,1);
permgeom        = master.permgeom(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);
shapgeomft      = master.shapmf(:,:,1);
dshapgeomft     = reshape(permute(master.shapmf(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);

%UH0 = permute(reshape(UH0,[ncu npf nf]),[2 1 3]);
elcon = reshape(mesh.elcon,[npf mesh.nfe mesh.ne]);
UH   = zeros(npf,ncu,nf);
bhg  = zeros(ngf,ncu);
for i = 1:nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face                
        kf = mesh.t2f(fi,:);         % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element
        i2 = kf(2,:)==i;  % obtain the index of face i in the 2nd element                                            
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        j2 = elcon(:,i2,fi(2)) - (i-1)*npf;                
                
        udg1 = UDG(perm(j1,i1),:,fi(1));
        udg2 = UDG(perm(j2,i2),:,fi(2));        
        udg1g = shapft*udg1;
        udg2g = shapft*udg2;
                
        p = mesh.dgnodes(permgeom(j1,i1),:,fi(1));             
        pg = shapgeomft*p;        
        dpg = dshapgeomft*p;                
        switch nd
            case 2
                jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
                nl  = [dpg(:,2)./jac,-dpg(:,1)./jac];                
            case 3
                dpg = permute(reshape(dpg,[npf nd-1 nd]), [1 3 2]);    
                nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
                nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
                nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
                jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
                nl  = bsxfun(@rdivide, nlg, jac);
            otherwise
                error('Dimension is not implemented');
        end                
        
        param = app.arg;
        time = app.time;        
        An1 = getbn(0,nl,udg1g,param);
        An2 = getbn(0,nl,udg2g,param);
        An  = An1+An2;
                
        H = reshape(shapfgdotshapfc*bsxfun(@times,An,jac),[npf npf ncu ncu]);        
        H = reshape(permute(H,[1 3 2 4]),[npf*ncu npf*ncu]);
        
        f1 = flux(pg,udg1g,param,time);
        f2 = flux(pg,udg2g,param,time);        
        df = f1-f2;        
        for  j = 1:ncu
            bhg(:,j) = sum(reshape(df(:,j,:),[ngf nd]).*nl,2) + ...
                       sum(reshape(An1(:,j,:),[ngf ncu]).*udg1g(:,1:ncu),2) + ...
                       sum(reshape(An2(:,j,:),[ngf ncu]).*udg2g(:,1:ncu),2);
        end
        uhg = shapfg*bsxfun(@times,bhg,jac);        
        
        UH(:,:,i) = H\uhg;         
        %[UH(:,:,i)-UH0(:,:,i)]
        %pause
    else % face i is a boundary face
        %UH(:,:,i) = U(perm(j1,i1),:,fi(1));        
        
        kf = mesh.t2f(fi(1),:); % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element                 
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        udg = UDG(perm(j1,i1),:,fi(1));        
        p = mesh.dgnodes(perm(j1,i1),:,fi(1));             
        b = -fi(2);
        ib = app.bcm(b);
        uinf = app.bcs(b,:);
        param = app.arg;
        time = app.time;                        
        
        % normal vector    
        pg = shapgeomft*p;        
        dpg = dshapgeomft*p;                
        switch nd
            case 2
                jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
                nl  = [dpg(:,2)./jac,-dpg(:,1)./jac];                
            case 3
                dpg = permute(reshape(dpg,[npf nd-1 nd]), [1 3 2]);    
                nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
                nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
                nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
                jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
                nl  = bsxfun(@rdivide, nlg, jac);
            otherwise
                error('Dimension is not implemented');
        end                
           
        udgg = shapft*udg;
        An = getbn(ib,nl,udgg,param);                
        H = reshape(shapfgdotshapfc*bsxfun(@times,An,jac),[npf npf ncu ncu]);        
        H = reshape(permute(H,[1 3 2 4]),[npf*ncu npf*ncu]);

        bhg = uhatbou(ib,uinf,nl,pg,udgg,param,time);                
        uhg = shapfg*bsxfun(@times,bhg,jac);   
        UH(:,:,i) = H\uhg;                     
    end        
end 
UH = reshape(permute(UH,[2 1 3]),[ncu npf*nf]);

%max(abs(UH(:)-UH0(:)))



