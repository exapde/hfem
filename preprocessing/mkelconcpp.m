function [elcon,facecon,edg] = mkelconcpp(p,t,f,t2f,elementtype,isEDGface,porder,dim,check)
%MKT2T Compute Element to Element Connectivity.
%   T2T=MKT2T(T,MESHTYPE)
%
%      T:         Triangle indices (NT,3)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%
%      T2T:       Triangle to Trangle Connectivity (NT,3)
%                 T2T(IT,IN) is the trangle that shares an edge
%                 with triangle IT and does nont contain node T(IT,IN).
%                 When an element is adjacent to the boundary the
%                 corresponding entry in T2T is set to zero
%
if nargin<=8
    check = 0;
end

elem = unique(elementtype);
nelem = length(elem);
nfe = zeros(nelem);
nve = nfe;
npe = nfe;
nvf = zeros(nelem,6);
npf = zeros(nelem,6);
phifc = []; 
phivl = [];
perm = [];
for i=1:nelem
    elemtype = elem(i);    
    [philocvl,philocfc,plocvl,plocfc,permface] = localbasis(porder,dim,elemtype); 
    phifc = [phifc; reshape(cell2mat(philocfc),[],1)];
    perm = [perm; reshape(cell2mat(permface),[],1)];
    phivl = [phivl; philocvl(:)];
    npe(i) = size(plocvl,1);
    for j=1:length(permface)
        npf(i,j) = length(permface{j});
    end
    switch dim
        case 1
            nfe(i) = 2;
            nve(i) = 2;
            nvf(i,1:2) = [1 1];
        case 2
            if elemtype==0
                nfe(i) = 3;
                nve(i) = 3;
                nvf(i,1:3) = [2 2 2];
            elseif elemtype==1
                nfe(i) = 4;
                nve(i) = 4;
                nvf(i,1:4) = [2 2 2 2];
            end
        case 3            
            if elemtype==0
                nfe(i) = 4;
                nve(i) = 4;
                nvf(i,1:4) = [3 3 3 3];
            elseif elemtype==1
                nfe(i) = 6;
                nve(i) = 8;
                nvf(i,1:6) = [4 4 4 4 4 4];
            elseif elemtype==2
                nfe(i) = 5;
                nve(i) = 6;
                nvf(i,1:5) = [3 3 4 4 4];
            elseif elemtype==3
                nfe(i) = 5;
                nve(i) = 5;
                nvf(i,:) = [4 3 3 3 3];
            end
        otherwise
            error('Only can handle dim=1, dim=2 or dim=3');
    end    
end

ndims = [dim,size(t,2),size(f,2),size(p,1),max(npe),max(nve),max(nfe),max(npf(:)),max(nvf(:)),nelem];
if check==1
    [elcon,facecon,edg] = mkelconc(p,t,f,t2f,elementtype,elem,isEDGface,phivl,phifc,perm-1,npe,nve,nfe,npf,nvf,ndims);
else
    [elcon,facecon] = mkelconc(p,t,f,t2f,elementtype,elem,isEDGface,phivl,phifc,perm-1,npe,nve,nfe,npf,nvf,ndims);
    edg = [];
end

if check == 1
    mesh.p = p;
    mesh.t = t'+1;
    mesh.f = f(2:end,:)'+1;
    mesh.f(:,end+1) = isEDGface;
    mesh.t2f = t2f'+1;
    mesh.nd = dim;
    mesh.porder = porder(1);
    mesh.elemtype = elem(1);
    mesh.perm = reshape(perm,[npf(1) nfe(1)]);
    [elcon2,~,edg2] = mkelcon(mesh);    
    if (max(abs(elcon2(:)-elcon(:)-1))) >1e-12
        error('elcon is incorrect');
    end
    if max(abs(edg2(:)-edg(:))) >1e-12
        error('edg is incorrect');
    end
end

