function [t2f,t2t,f,ne,nf,nfe,nve,nvf] = mkt2fmatlab(t,elementtype,dim)

elem = unique(elementtype);
nelem = length(elem);

if (length(nelem)==1) && ismember(elem(1),[0 1])
    %[t2f,t2t,f,ne,nf,nfe,nve,nvf] = mkt2f(t,elem(1));    
    [t2f,t2t,f,ne,nf,nfe,nve,nvf] = mkt2f2(t,elem(1));    
    return;
end

nfe = zeros(nelem,1);
nve = nfe;
nvf = nfe;
nes = nfe;
face = cell(nelem,1);
nvfs = cell(nelem,1);
for i=1:nelem
    elemtype = elem(i);    
    nes(i) = length(find(elementtype==elemtype));    
    switch dim
        case 1
            nfe(i) = 2;
            nve(i) = 2;
            nvf(i) = 1;
            face{i} = [1;2];
            nvfs{i} = [1 1];
        case 2
            if elemtype==0
                nfe(i) = 3;
                nve(i) = 3;
                nvf(i) = 2;
                face{i} = [[2,3];[3,1];[1,2]];
                nvfs{i} = [2 2 2];
            elseif elemtype==1
                nfe(i) = 4;
                nve(i) = 4;
                nvf(i) = 2;
                face{i} = [[1,2];[2,3];[3,4];[4,1]];
                nvfs{i} = [2 2 2 2];
            end
        case 3            
            if elemtype==0
                nfe(i) = 4;
                nve(i) = 4;
                nvf(i) = 3;
                face{i} = [[2,3,4];[1,4,3];[1,2,4];[1,3,2]];
                nvfs{i} = [3 3 3 3];
            elseif elemtype==1
                nfe(i) = 6;
                nve(i) = 8;
                nvf(i) = 4;
                face{i} = [[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]];
                nvfs{i} = [4 4 4 4 4 4];
            elseif elemtype==2
                nfe(i) = 5;
                nve(i) = 6;
                nvf(i) = 4;                
                face{i} = [[1,3,2,0];[4,5,6,0];[2,3,6,5];[3,1,4,6];[1,2,5,4]];
                nvfs{i} = [3 3 4 4 4];
            elseif elemtype==3
                nfe(i) = 5;
                nve(i) = 5;
                nvf(i) = 4;                           
                face{i} = [[1,4,3,2];[1,2,5,0];[2,3,5,0];[3,4,5,0];[4,1,5,0]];  
                nvfs{i} = [4 3 3 3 3];
            end
        otherwise
            error('Only can handle dim=1, dim=2 or dim=3');
    end    
end

ne = size(t,1);
ns = sum(nfe.*nes);
nfemax = max(nfe);
nvfmax = max(nvf);
faces = zeros(ns,nvfmax);
facei = zeros(ns,2);
jf = 1;
for ie=1:ne
    elemtype = elementtype(ie);    
    i = find(elem==elemtype);
    nfei = nfe(i);    
    for j = 1:nfei
        nvfi = nvfs{i}(j);
        faces(jf,1:nvfi) = t(ie,face{i}(j,1:nvfi));
        facei(jf,:) = [ie j];        
        jf = jf+1;
    end
end

[f,~,jx]=unique(sort(faces,2,'descend'),'rows'); 
nf = size(f,1);
f = [0*f(:,1) f 0*f(:,1) 0*f(:,1)];
t2f = zeros(ne,nfemax);
t2t = zeros(ne,nfemax);
for i = 1:nf
    % find k to match faces to f
    k = find(jx==i);
    % convert k into element index e and local face index l
    e = facei(k,1);
    l = facei(k,2);    
    % compute f and t2f
    for m=1:length(k)
        f(i,nvfmax+m) = e(m);                
        t2f(e(m),l(m)) = i;        
    end
    nvfi = nvfs{e(1)}(l(1));
    f(i,1) = nvfi;
    f(i,2:nvfmax+1) = faces(k(1),:);
    % compute t2t
    if m>1
        t2t(e(1),l(1)) = e(2);
        t2t(e(2),l(2)) = e(1);
    end        
end


function [t2f,t2t,f,ne,nf,nfe,nve,nvf] = mkt2f(t,elemtype)
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

% npv : number of nodes per volume element
% nfv : number of faces per volume element
% npf : number of nodes per face element

if nargin<2, elemtype=0; end

[nt,npv]=size(t);

if npv==2
    dim=1;
    nfv=2;
else
    if elemtype==0 % tri/tet elements
        dim=size(t,2)-1;    
        nfv=dim+1;        
    elseif elemtype==1 % quad/hex elements
        dim=log2(size(t,2));   
        nfv=2*dim;
    end
end

switch dim
    case 1
        faces=[t(:,1)
               t(:,2)];
         nve = 2;  
    case 2
        if elemtype==0
            faces=[t(:,[2,3])
                   t(:,[3,1])
                   t(:,[1,2])];
            nve = 3;  
        elseif elemtype==1
            faces=[t(:,[1,2])
                   t(:,[2,3])
                   t(:,[3,4])
                   t(:,[4,1])];
            nve = 4;   
        end
    case 3
        if elemtype==0
            faces=[t(:,[2,3,4])
                   t(:,[1,4,3])
                   t(:,[1,2,4])
                   t(:,[1,3,2])];
             nve = 4;  
        elseif elemtype==1
            faces=[t(:,[1,4,3,2])
                   t(:,[5,6,7,8])
                   t(:,[1,2,6,5])
                   t(:,[3,4,8,7])
                   t(:,[2,3,7,6])
                   t(:,[4,1,5,8])];
            nve = 8;   
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end

nvf = size(faces,2);
nfe = nfv;
ne  = nt;
faces = reshape(faces,[ne nfe nvf]);
faces = reshape(permute(faces,[2 1 3]),[nfe*ne nvf]);

[f,~,jx]=unique(sort(faces,2),'rows'); 
nf = size(f,1);
f = [f 0*f(:,1) 0*f(:,1)];
t2f = zeros(ne,nfe);
t2t = zeros(ne,nfe);
for i = 1:nf
    % find k to match faces to f
    k = find(jx==i);
    % convert k into element index e and local face index l
    e = ceil(k/nfe);
    l = k - (e-1)*nfe;    
    % compute t2f
    for m=1:length(k)
        f(i,nvf+m) = e(m);                
        t2f(e(m),l(m)) = i;        
    end
    f(i,1:nvf) = faces(k(1),:);
    % compute t2t
    if m>1
        t2t(e(1),l(1)) = e(2);
        t2t(e(2),l(2)) = e(1);
    end        
end
f = [nvf*ones(nf,1) f];


function [t2f,t2t,f,ne,nf,nfe,nve,nvf] = mkt2f2(t,elemtype)
%MKT2F Compute Face Connectivity and Triangle to Face Connetivity.
%   [F,T2F]=MKT2F(T)
%
%      T:         Triangle indices (NT,3)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%
%      F:         Face connectivity (NF,4) (for boundary edges F(:,4)=0)
%      T2F:       Triangle to Face Connectivity (NT,3)
%
%   See also: MKT2T.

% npv : number of nodes per volume element
% nfv : number of faces per volume element
% npf : number of nodes per face element

if nargin<2, elemtype=0; end

[nt,npv]=size(t);
if npv==2 % 1D
    dim=1;
    nfv=2;
    npf=1;
else
    if elemtype==0 % tri/tet elements
        dim=size(t,2)-1;        
        nfv=dim+1;
        npf = dim;
    else % quad/hex elements
        dim=log2(size(t,2));        
        nfv=2*dim;
        npf=2*(dim-1);
    end
end

switch dim
    case 1
        face=[1;2];
        nve = 2;   
        nvf = 1;
    case 2
        if elemtype==0
            face=[[2,3];[3,1];[1,2]];
            nve = 3;
            nvf = 2;
        elseif elemtype==1
            face=[[1,2];[2,3];[3,4];[4,1]];
            nve = 4;
            nvf = 2;
        end
    case 3
        if elemtype==0
            face=[[2,3,4];[1,4,3];[1,2,4];[1,3,2]];
            nve = 4;
            nvf = 3;
        elseif elemtype==1
            face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]];
            nve = 8;   
            nvf = 4;
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end


t2t = mkt2t(t,elemtype);
nb = sum(sum(t2t <= 0));
f = zeros((nfv*nt+nb)/2,npf+2);
t2f = zeros(nt,nfv);
jf = 0;
for i=1:nt
    for j=1:nfv
        if t2t(i,j) > i || t2t(i,j) <=0
            ie = t2t(i,j);
            jf = jf + 1;
            
            f(jf,1:npf) = t(i,face(j,:));
            f(jf,npf+1) = i;
            f(jf,npf+2) = ie;
            t2f(i,j) = jf;
            
            if ie > 0
                %k = sum(reshape(t(ie,face),[nfv npf]),2)-sum(f(jf,1:npf))==0;
                %t2f(ie,k) = jf;
                a = sort(reshape(t(ie,face),[nfv npf]),2);
                b = sort(f(jf,1:npf));
                k = sum(abs(a-repmat(b,[nfv 1])),2)==0;
                t2f(ie,k) = jf;
            end                        
        end
    end
end

nfe = nfv;
ne = nt;
nf = size(f,1);
f = [nvf*ones(nf,1) f];





