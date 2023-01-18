
function [S_j, I_j, Itilde_j] = jumpShockSensor(mesh,master,UDG,hFlag,fieldNo,smoothSensorFlag, plotFlag)

if nargin < 4; hFlag = 1; end
if nargin < 5; fieldNo = 1; end
if nargin < 6; smoothSensorFlag = 0; end
if nargin < 7; plotFlag = 0; end

if hFlag == 2; warning('Using "h = minimum element height" in the face jump-based shock sensor is not recommended.'); end

I_j_lim = 1;
k = 50;

nd = mesh.nd;
ne = mesh.ne;
nf = mesh.nf;
npv = master.npv;
npf = master.npf;
ngv = master.ngv;
porder = mesh.porder;
perm = mesh.perm;
elcon = mesh.elcon;
dgnodes = mesh.dgnodes;

shapft = master.shapft(:,:,1);
shapfg = master.shapfg(:,:,1);
shapmf = master.shapmf;

shapvt = squeeze(master.shapvt(:,:,1));
shapmv = master.shapmv;
gwvl = master.gwvl;

absJump = zeros(ne,1);
faceMeasure = zeros(ne,1);
mapping = zeros(npf,1);

[~, ~, jacE] = volgeom(shapmv,permute(dgnodes(:,1:nd,:),[1,3,2]));
jacE = reshape(jacE, [ngv,ne]);

for face = 1:nf
    if mesh.f(face,end) > 0     % Interior face
        elem1 = mesh.f(face,end-1);
        elem2 = mesh.f(face,end);
        
        localFace1 = find(mesh.t2f(elem1,:) == face);
        if length(localFace1) ~= 1; error('Something wrong.'); end
        localFace2 = find(mesh.t2f(elem2,:) == face);
        if length(localFace2) ~= 1; error('Something wrong.'); end
        
        UHnodes1 = elcon((localFace1-1)*npf+1:localFace1*npf,elem1);
        UHnodes2 = elcon((localFace2-1)*npf+1:localFace2*npf,elem2);
        for i = 1:npf
            tmp = find(UHnodes2 == UHnodes1(i));
            if length(tmp) ~= 1; error('Something wrong.'); end
            mapping(i) = tmp;
        end
        
        fields1 = shapft*UDG(perm(:,localFace1),:,elem1);                
        fields2 = shapft*UDG(perm(mapping,localFace2),:,elem2); 
        
        field1 = getField(fields1,fieldNo);
        field2 = getField(fields2,fieldNo);
        
        xDG1 = mesh.dgnodes(perm(:,localFace1),1:nd,elem1);
        
        [~,~,jacF] = facegeom(shapmf,xDG1);
        
        intNum = sum(shapfg*diag(jacF)*abs(field1-field2));
        intDen = sum(shapfg*jacF(:));
        
        absJump(elem1) = absJump(elem1) + intNum;
        absJump(elem2) = absJump(elem2) + intNum;
        
        faceMeasure(elem1) = faceMeasure(elem1) + intDen;
        faceMeasure(elem2) = faceMeasure(elem2) + intDen;
    else
        % Boundary face
    end
end

udggField = shapvt*squeeze(getField(UDG,fieldNo));

% mesh = getElementSize(mesh);
% hField_g = shapvt*squeeze(mesh.hField);

fieldAvg = zeros(ne,1);
% hElem = zeros(ne,1);
for elem=1:ne
    fieldAvg(elem) = sum(gwvl .* udggField(:,elem) .* jacE(:,elem)) / sum(gwvl .*jacE(:,elem));
%     hElem(elem) = sum(gwvl .* hField_g(:,elem) .* jacE(:,elem)) / sum(gwvl .*jacE(:,elem));
end
fieldAvg = fieldAvg(:);
% hElem = hElem(:);

if hFlag == 0       % Use constant h = 0.02
    hElem = 0.02 * ones(ne,1);
elseif hFlag == 1   % Use average h = measure^(1/nd)
    measure = computeElemMeasure(mesh,master);
    hElem = measure.^(1/nd);
    hElem = hElem(:);
elseif hFlag == 2   % Use h = minimum height
    mesh = getElementSize(mesh);
    hField_g = shapvt*squeeze(mesh.hField);

    hElem = zeros(ne,1);
    for elem=1:ne
        hElem(elem) = sum(gwvl .* hField_g(:,elem) .* jacE(:,elem)) / sum(gwvl .*jacE(:,elem));
    end
    hElem = hElem(:);
end

I_j = absJump ./ (faceMeasure .* hElem.^((porder+1)/2) .* fieldAvg);
Itilde_j = absJump ./ (faceMeasure .* fieldAvg);

I_j_p1CG = p0DG_to_p1CG(I_j, mesh);

I_j = repmat(I_j(:)',[npv,1]);
Itilde_j = repmat(Itilde_j(:)',[npv,1]);

if smoothSensorFlag == 0
    S_j = 1 ./ (1 + exp( - k *(I_j-I_j_lim)));
elseif smoothSensorFlag == 1
    S_j = 1 ./ (1 + exp( - k *(I_j_p1CG-I_j_lim)));
end

if plotFlag
    figure(1); scaplot(mesh,S_j); figure(2); scaplot(mesh,log(I_j));
    figure(3); scaplot(mesh,log(I_j_p1CG)); figure(4); scaplot(mesh,log(Itilde_j));
end

end


function field = getField(fields,fieldNo)

if isa(fieldNo,'numeric')
    if floor(fieldNo) ~= ceil(fieldNo); error('Something wrong.');
    else
        field = fields(:,fieldNo,:);
    end
elseif strcmp(fieldNo,'s')
    error('Type of field not implemented.');
else
    error('Type of field not implemented.');
end

end


function [pg, nlg, jac] = facegeom(shapgeomft,pn)
% FACEGEOM computes dg nodes, Jacobian determinant and normal vectors at Gauss points  

%   [pg, nlg, jac] = facegeom(shapgeomft,dgnodes,perm)
%
%    SHAPGEOMFT :  Shape functions and derivatives at Gauss points
%    DGNODES    :  Geometry DG nodes 
%    PERM       :  Indices of the boundary nodes
%    PG         :  Physical nodes at Gauss points 
%    nlg        :  Normal vector at Gauss points
%    jac        :  Determinant of the Jacobian mapping 

nq    = size(pn,2);
ngf   = size(shapgeomft,1);
npf   = size(shapgeomft,2);
nd    = size(shapgeomft,3);

if nd>1
    dshapft  = reshape(permute(shapgeomft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
    dpg = dshapft*pn(:,1:nd);
    dpg = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 4 5 2]);    
    dpg = reshape(dpg,[ngf,nd,nd-1]);    
end

shapgeomft   = shapgeomft(:,:,1);
pg = shapgeomft*pn;
pg = reshape(pg,[ngf nq]);

switch nd
    case 1
        jac = ones(1,1);
        nlg = [ones(1,1); -ones(1,1)];        
        nlg = nlg(:);
    case 2
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
    case 3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
    otherwise
        error('Dimension is not implemented');
end

end
