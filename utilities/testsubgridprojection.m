function testsubgridprojection

pgeom = 4;
ppr   = 4;
pgausspr = 8;
nrefpr = 2;
ndim = 2;
elemtype = 0;
nodetype = 0;

pdu   = 2;
pgaussdu = 4;
nrefdu = 2;

masterpr = mkmastersubgrid(pgeom,ppr,pgausspr,nrefpr,ndim,elemtype,nodetype);
masterdu = mkmastersubgrid(pgeom,pdu,pgaussdu,nrefdu,ndim,elemtype,nodetype);

ne = 10;
nc = 4;
npm = size(masterpr.shapmv,2);
ntpr = size(masterpr.t,1);
npvpr = size(masterpr.shapvl,1);

dgnodes = rand(npm,ndim,ne);
udgpr = rand(npvpr,nc,ntpr,ne);

[udgdu, dgmshapdu, dgmshappr, pdmap, gaussprdu, shapprdu] = subgridprojection(masterpr,masterdu,dgnodes,udgpr);

gridtype = [4 0; 2 1; 0 2];
mesh.gridtype = gridtype;
mesh.subgrids = randi(size(gridtype,1),ne,1);
mesh.dgnodes = dgnodes;
mesh.nd = ndim;
mesh.nodetype = nodetype;
mesh.elemtype = elemtype;
mesh.porder = gridtype(1,1);

for i = 1:size(mesh.gridtype,1)
    mastersubgrid{i} = mkmastersubgrid(mesh.porder,mesh.gridtype(i,1),max(2*mesh.gridtype(i,1),1),mesh.gridtype(i,2),mesh.nd,mesh.elemtype,mesh.nodetype);        
    npv  = size(mastersubgrid{i}.plocvl,1);
    npf  = size(mastersubgrid{i}.perm,1);        
    nf   = max(mastersubgrid{i}.elcon(:))/(npf);
    nes  = size(mastersubgrid{i}.t,1);    
    ind  = find(mesh.subgrids==i);
    nt   = length(ind);
    udg{i} = rand(npv,nc,nes*nt);
end

upg = udgproj(mesh,mastersubgrid,udg);

newsg = randi(size(gridtype,1),ne,1);
[udgt,uhat] = subgridadaptation(mesh,mastersubgrid,udg,newsg,nc);

save projectiontest.mat udgdu udgpr dgnodes masterpr masterdu ...
     pgeom ppr pgausspr nrefpr ndim elemtype nodetype pdu pgaussdu nrefdu ...
     dgmshapdu dgmshappr pdmap gaussprdu shapprdu mesh mastersubgrid udg upg udgt uhat newsg
 
 
function [udgt,uhat] = subgridadaptation(mesh,mastersubgrid,udg,newsg,nc)

ns = length(mastersubgrid);

oldsg = mesh.subgrids;
ds = oldsg - newsg;

%es = find(ds~=0); % list of modified elements    
for i=1:ns
    indi{i} = find(newsg==i);
    indj{i} = find(oldsg==i);
end    

for i=1:ns % compute projection on subgrid i
    if isempty(indi{i})==0
        ni  = length(indi{i});
        ngi = size(mastersubgrid{i}.elemnodes,3); 
        npi = size(mastersubgrid{i}.shapvl,1); 
        udgt{i} = zeros(npi,nc,ngi*ni);
        for j=1:ns
            s{j} = [];
        end
        for k=1:ni
            e = indi{i}(k); % element e                           
            if ds(e)==0
                m = find(indj{i}==e);
                udgt{i}(:,:,(k-1)*ngi+1:k*ngi) = udg{i}(:,:,(m-1)*ngi+1:m*ngi);                                        
            else
                j   = oldsg(e);          % subgrid j                  
                m   = find(indj{j}==e);  % position of element e in the list sj                                        
                s{j} = [s{j}; e k m];                        
            end
        end
        for j=1:ns
            if isempty(s{j})==0                    
                ngj = size(mastersubgrid{j}.elemnodes,3); % number of sub-elements for subgrid type j                                
                npj = size(mastersubgrid{j}.shapvl,1); 
                nj1 = size(s{j},1);
                k1  = repmat((s{j}(:,2)'-1)*ngi,[ngi 1]) + repmat((1:ngi)',[1 nj1]);                
                k2  = repmat((s{j}(:,3)'-1)*ngj,[ngj 1]) + repmat((1:ngj)',[1 nj1]);
                tm  = subgridprojection(mastersubgrid{j},mastersubgrid{i},mesh.dgnodes(:,:,s{j}(:,1)),reshape(udg{j}(:,:,k2),[npj nc ngj nj1]));                                                                                                                                
                udgt{i}(:,:,k1) = reshape(tm,[npi nc ngi*nj1]);                
            end
        end
    end
end    


ncu = nc;
for i=1:ns
    npf  = size(mastersubgrid{i}.perm,1);        
    nf   = max(mastersubgrid{i}.elcon(:))/(npf);
    nes  = size(mastersubgrid{i}.t,1);    
    ind  = find(newsg==i);
    ne   = length(ind);    
    uhat{i} = zeros(ncu,npf*nf,ne);           
    for k = 1:ne                
        uhat{i}(:,:,k) = inituhat(mastersubgrid{i},mastersubgrid{i}.elcon,udgt{i}(:,:,(k-1)*nes+1:k*nes),ncu);  
    end    
end
 
function udg = udgproj(mesh,mastersubgrid,u)

ns = max(mesh.subgrids);
npm = size(mesh.dgnodes,1);
ne  = size(mesh.dgnodes,3);

for j=1:ns
    nc  = size(u{j},2);
    if nc > 0
        break;
    end
end

udg = zeros(npm,nc,ne);
%mesht = mesh;
for j = 1:ns     
    ind = find(mesh.subgrids==j);    
    if isempty(ind)==0
        if j==1
            udg(:,:,ind) = u{j};
        else
            nj = length(ind);
            ng = size(mastersubgrid{j}.elemnodes,3);        
            npv = size(u{j},1);
            %mesht.dgnodes = mesh.dgnodes(:,:,ind);                                             
            udg(:,:,ind) = subgridprojection(mastersubgrid{j},mastersubgrid{1},mesh.dgnodes(:,:,ind),reshape(u{j},[npv nc ng nj]));                            
        end
    end
end
 
function [udgdu, dgmshapdu, dgmshappr, pdmap, gaussprdu, shapprdu] = subgridprojection(mastergridpr,mastergriddu,dgnodes,udgpr)
%SUBGRIDPROJECTION perform the L2 projection from the primal space to the
%dual space
%
%   udgdu = SUBGRIDPROJECTION(mastergridpr,mastergriddu,dgnodes,udgpr)
%
%    MASTERGRIDPR: Mastergrid structure for the primal space
%    MASTERGRIDDU: Mastergrid structure for the dual space
%    DGNODES:      Geometry DG nodes 
%    UDGPR:        Field variables on the primal space
%    UDGDU:        Field variables on the dual space

ne = size(dgnodes,3);
nd = size(dgnodes,2);
nc = size(udgpr,2);
n = mastergridpr.nref - mastergriddu.nref;
ntdu     = mastergriddu.nt;      % number of elements in the dual subgrid
ntpr     = mastergridpr.nt;      % number of elements in the primal subgrid
elemtype = mastergridpr.elemtype;
nodetype = mastergridpr.nodetype;

% porderpr = mastergridpr.porder; % order for the primal subgrid
% nrefpr   = mastegridpr.nref;    % refinement level for the primal subgrid
% ntpr     = mastegridpr.nt;      % number of elements in the primal subgrid
% plocvlpr = mastergridpr.plocvlpr;
% elemtype = mastergridpr.elemtype;
% nodetype = mastergridpr.nodetype;
% 
% porderdu = mastergriddu.porder; % order for the dual subgrid
% nrefdu   = mastegriddu.nref;    % refinement level for the dual subgrid
% ntdu     = mastegriddu.nt;      % number of elements in the dual subgrid
% npvdu    = mastegriddu.npv;     % number of basis functions per element in the dual subgrid
% shapgeomdu = mastergriddu.shapgeom;
% shapvldu  = mastergriddu.shapvl;
% gpvldu  = mastergriddu.gpvl;
% gwvldu  = mastergriddu.gwvl;

% derivatives of geometry shape functions on the dual element
dgmshapdu = zeros(mastergriddu.ngv,mastergriddu.npm,nd);
for i=2:nd+1
    dgmshapdu(:,:,i-1) = mastergriddu.shapmv(:,:,i);
end
dgmshapdu = reshape(permute(dgmshapdu,[1 3 2]),[mastergriddu.ngv*nd mastergriddu.npm]);

% derivatives of geometry shape functions on the primal element
dgmshappr = zeros(mastergridpr.ngv,mastergridpr.npm,nd);
for i=2:nd+1
    dgmshappr(:,:,i-1) = mastergridpr.shapmv(:,:,i);
end
dgmshappr = reshape(permute(dgmshappr,[1 3 2]),[mastergridpr.ngv*nd mastergridpr.npm]);

udgdu = zeros(mastergriddu.npv,nc,ntdu,ne);
if (n==0) % the primal and dual subgrids are the same   
    pdmap=0;
    gaussprdu=0;
    
    % primal shape functions at the Gauss points of the dual element
    shapprdu = mkshape(mastergridpr.porder,mastergridpr.plocvl,mastergriddu.gpvl,elemtype);
    
    for i=1:ne % for each subgrid
        for k=1:ntdu % for each element on the dual subgrid
            % compute subgrid DG nodes in the physical space    
            dgx = mastergriddu.shapgeom(:,:,k)*dgnodes(:,:,i);                           

            % compute the Jacobian matrix at Gauss points
            Jg = dgmshapdu*dgx(:,1:nd);
            Jg = reshape(Jg,[mastergriddu.ngv nd nd]);        
            jac = volgeom(Jg);   
            
            % mass matrix on the dual element
            M = mastergriddu.shapvl(:,:,1)*diag(mastergriddu.gwvl.*jac)*(mastergriddu.shapvl(:,:,1))';            
            
            % cross mass matrix
            N = mastergriddu.shapvl(:,:,1)*diag(mastergriddu.gwvl.*jac)*(shapprdu(:,:,1))';            
            
            % projection
            udgdu(:,:,k,i) = M\(N*udgpr(:,:,k,i));            
        end        
    end
elseif n>0 % the primal subgrid has more elements than dual subgrid
    
    % primal-dual grid
    [ps,ts] = masternodes(max(2^abs(n),1),nd,elemtype,nodetype);
    nts = size(ts,1);
    
    pdmap = zeros(ntdu,nts);
    for i=1:ntdu % for each dual element 
        for j=1:ntpr % for each primal element 
            % obtain the vertex points for the primal element
            pr = mastergridpr.p(mastergridpr.t(j,:),:);
            for k=1:nts
                pd = mapp(ps(ts(k,:),:),mastergriddu.p(mastergriddu.t(i,:),:));                    
                if norm(mean(pr)-mean(pd)) <= 1e-6
                    pdmap(i,k) = j;                    
                end
            end
        end
    end    
    
    % gauss nodes for the primal-dual subgrid
    gaussprdu = mknodes(ps,ts,mastergridpr.gpvl);
    
    % dual shape functions at the Gauss points of the primal-dual subgrid    
    for j=1:nts
        tmp = mkshape(mastergriddu.porder,mastergriddu.plocvl,gaussprdu(:,:,j),elemtype);
        shapprdu(:,:,j) = tmp(:,:,1);
    end
    
    for i=1:ne % for each subgrid
        for k=1:ntdu % for each element on the dual subgrid
            % compute subgrid DG nodes in the dual subgrid
            dgx = mastergriddu.shapgeom(:,:,k)*dgnodes(:,:,i);                           

            % compute the Jacobian matrix at Gauss points
            Jg = dgmshapdu*dgx(:,1:nd);
            Jg = reshape(Jg,[mastergriddu.ngv nd nd]);        
            jac = volgeom(Jg);   
            
            % mass matrix on the dual element
            M = mastergriddu.shapvl(:,:,1)*diag(mastergriddu.gwvl.*jac)*(mastergriddu.shapvl(:,:,1))';            
            
            for j=1:nts       
                % element index of the primal subgrid
                a = pdmap(k,j);
                
                % compute subgrid DG nodes in the primal subgrid
                dgpr = mastergridpr.shapgeom(:,:,a)*dgnodes(:,:,i);                           

                % compute the Jacobian matrix at Gauss points
                Jgpr = dgmshappr*dgpr(:,1:nd);
                Jgpr = reshape(Jgpr,[mastergridpr.ngv nd nd]);        
                jacpr = volgeom(Jgpr);   
                
                % cross mass matrix                
                N = shapprdu(:,:,j)*diag(mastergridpr.gwvl.*jacpr)*(mastergridpr.shapvl(:,:,1))';            
                                
                if j==1
                    F = N*udgpr(:,:,a,i);
                else
                    F = F + N*udgpr(:,:,a,i);
                end
            end
            
            % projection
            udgdu(:,:,k,i) = M\F;            
        end        
    end    
else % the primal subgrid has less elements than dual subgrid
    % primal-dual grid
    [ps,ts] = masternodes(max(2^abs(n),1),nd,elemtype,nodetype);
    nts = size(ts,1);
    
    pdmap = zeros(ntdu,2);
    for i=1:ntdu % for each dual element 
        % obtain the vertex points for the dual element
        pd = mastergriddu.p(mastergriddu.t(i,:),:);
        for j=1:ntpr % for each primal element             
            for k=1:nts
                pr = mapp(ps(ts(k,:),:),mastergridpr.p(mastergridpr.t(j,:),:));                    
                if norm(mean(pr)-mean(pd)) <= 1e-6
                    pdmap(i,1) = j;                    
                    pdmap(i,2) = k;                    
                end
            end
        end
    end            
    
    % gauss nodes for the primal-dual subgrid
    gaussprdu = mknodes(ps,ts,mastergriddu.gpvl);
    
    % primal shape functions at the Gauss points of the primal-dual subgrid    
    for j=1:nts
        tmp = mkshape(mastergridpr.porder,mastergridpr.plocvl,gaussprdu(:,:,j),elemtype);
        shapprdu(:,:,j) = tmp(:,:,1);
    end   
    
    for i=1:ne % for each subgrid
        for k=1:ntdu % for each element on the dual subgrid
            % compute subgrid DG nodes in the dual subgrid
            dgx = mastergriddu.shapgeom(:,:,k)*dgnodes(:,:,i);                           

            % compute the Jacobian matrix at Gauss points
            Jg = dgmshapdu*dgx(:,1:nd);            
            Jg = reshape(Jg,[mastergriddu.ngv nd nd]);        
            jac = volgeom(Jg);               
            
            % mass matrix on the dual element
            M = mastergriddu.shapvl(:,:,1)*diag(mastergriddu.gwvl.*jac)*(mastergriddu.shapvl(:,:,1))';            
            
            % element index of the primal subgrid
            a = pdmap(k,1);
            b = pdmap(k,2);

            % cross mass matrix
            N = mastergriddu.shapvl(:,:,1)*diag(mastergriddu.gwvl.*jac)*(shapprdu(:,:,b))';                                                            
                
            % projection
            udgdu(:,:,k,i) = M\(N*udgpr(:,:,a,i));            
        end        
    end        
end


function [jac,Xx] = volgeom(Jg)

ngv = size(Jg,1); 
nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;
        Xx = -ones(ngv,1);
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);
        Xx(:,1,1) = -Jg(:,2,2);
        Xx(:,2,1) = Jg(:,2,1);
        Xx(:,1,2) = Jg(:,1,2);
        Xx(:,2,2) = -Jg(:,1,1);
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);            
        Xx(:,1,1) = Jg(:,2,3).*Jg(:,3,2) - Jg(:,2,2).*Jg(:,3,3);
        Xx(:,2,1) = Jg(:,2,1).*Jg(:,3,3) - Jg(:,2,3).*Jg(:,3,1);
        Xx(:,3,1) = Jg(:,2,2).*Jg(:,3,1) - Jg(:,2,1).*Jg(:,3,2);
        Xx(:,1,2) = Jg(:,1,2).*Jg(:,3,3) - Jg(:,1,3).*Jg(:,3,2);
        Xx(:,2,2) = Jg(:,1,3).*Jg(:,3,1) - Jg(:,1,1).*Jg(:,3,3);
        Xx(:,3,2) = Jg(:,1,1).*Jg(:,3,2) - Jg(:,1,2).*Jg(:,3,1);
        Xx(:,1,3) = Jg(:,1,3).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,3);
        Xx(:,2,3) = Jg(:,1,1).*Jg(:,2,3) - Jg(:,1,3).*Jg(:,2,1);
        Xx(:,3,3) = Jg(:,1,2).*Jg(:,2,1) - Jg(:,1,1).*Jg(:,2,2);
    otherwise
        error('Dimension is not implemented');
end


function p=mapp(p,x)

nx = size(x,1);
nd = size(x,2);

if nd==1
     r=[0,1];
     
     C=x(:,1)'/[1,1; r];
     
     p=[0*p+1,p]*C';
elseif (nd==2) && (nx==4) % 2-D bilinear map    
     r=[0,1,1,0];
     s=[0,0,1,1];
     
     C=x(:,1)'/[1,1,1,1; r; r.*s; s];
     D=x(:,2)'/[1,1,1,1; r; r.*s; s];

     px=p(:,1); py=p(:,2);
     p=[0*px+1, px, px.*py, py]*[C;D]';
elseif (nd==3) && (nx==8) % 3-D bilinear map    
     r=[0,1,1,0,0,1,1,0];
     s=[0,0,1,1,0,0,1,1];
     t=[0,0,0,0,1,1,1,1];
     
     C=x(:,1)'/[1,1,1,1,1,1,1,1; r; r.*s; s; t; r.*t; r.*s.*t; s.*t];
     D=x(:,2)'/[1,1,1,1,1,1,1,1; r; r.*s; s; t; r.*t; r.*s.*t; s.*t];
     E=x(:,3)'/[1,1,1,1,1,1,1,1; r; r.*s; s; t; r.*t; r.*s.*t; s.*t];
             
     px=p(:,1); py=p(:,2); pz=p(:,3);
     p=[0*px+1, px, px.*py, py, pz, px.*pz, px.*py.*pz, py.*pz]*[C;D;E]';
elseif (nd==2) && (nx==3) % 2-D affine map    
     r=[0,1,0];
     s=[0,0,1];
     
     C=x(:,1)'/[1,1,1; r; s];
     D=x(:,2)'/[1,1,1; r; s];
     px=p(:,1); py=p(:,2);
     p=[0*px+1, px, py]*[C;D]'; 
elseif (nd==3) && (nx==4) % 3-D affine map    
     r=[0,1,0,0];
     s=[0,0,1,0];
     t=[0,0,0,1];
     
     C=x(:,1)'/[1,1,1,1; r; s; t];
     D=x(:,2)'/[1,1,1,1; r; s; t];
     E=x(:,3)'/[1,1,1,1; r; s; t];

     px=p(:,1); py=p(:,2); pz=p(:,3);
     p=[0*px+1, px, py, pz]*[C;D;E]';
end















