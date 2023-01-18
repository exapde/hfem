function UH = ldg_uhat(mesh,master,app,U)

ncu = size(U,2);   % number of components of UH
nf  = size(mesh.f,1);    % number of faces
npf = size(mesh.perm,1); % number of points per face

ubou = str2func(app.ubou);
shapfc = mkshape(mesh.porder,master.plocfc,master.plocfc,mesh.elemtype);

elcon = reshape(mesh.elcon,[npf mesh.nfe mesh.ne]);
perm = mesh.perm;
UH   = zeros(npf,ncu,nf);
for i = 1:nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face                
        kf = mesh.t2f(fi,:);         % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element
        i2 = kf(2,:)==i;  % obtain the index of face i in the 2nd element                                            
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        j2 = elcon(:,i2,fi(2)) - (i-1)*npf;                
        UH(:,:,i) = 0.5*(U(perm(j1,i1),:,fi(1)) + U(perm(j2,i2),:,fi(2)));                        
    else % face i is a boundary face
        %UH(:,:,i) = U(perm(j1,i1),:,fi(1));        
        
        kf = mesh.t2f(fi(1),:); % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element                 
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        u = U(perm(j1,i1),:,fi(1));        
        p = mesh.dgnodes(perm(j1,i1),:,fi(1));             
        b = -fi(2);
        ib = app.bcm(b);
        uinf = app.bcs(b,:);
        param = app.arg;
        time = app.time;
        
        % normal vector
        nd = mesh.nd;
        dshapft  = reshape(permute(shapfc(:,:,2:nd),[1 3 2]),[npf*(nd-1) npf]);
        dpg = dshapft*p;                
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
                
        UH(:,:,i) = ubou(ib,uinf,nl,p,u,param,time);                
    end        
end 

