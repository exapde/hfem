
% TODO:
% - Compare actual N to N given by solving
% Orr-Sommerfeld equation. (Does Drela have a solver?)
% - Compute stagnation point location and use it whenever necessary.

% DOUBTS:
% 1. ue = u*(ne) or u(ne)? (Alejandra does u*(ne))
% Indeed, she ALWAYS uses the pseudo-velicity u*.
% 3. How are A10 and A20 (amplifications at onset of transition) defined?

function [] = LESpostProcess(fileName, numProc, ib, mesh, master, app, dmd)

% ib: Airfoil boundary
% fileName = '3DLESIEDGp2h2';
% numPro: Number of processors used to compute numerical solution

% load([fileName, '.mat']);

% alpha = 0;
alpha = -0.527999463;       % T106C
% alpha = 0.493404579538797;      % Angle between aifoil chord and x-axis
LEcoord = [0, 0];               % (x,y) coordinates of leading edge

timeStepStart = 1100;
timeStepEnd = 2100;
timeStepFreq = 2;

integralAccuracy = 1;

% Parameters for envelope method transition prediction:
N1_cr_lowerSide = 6.0;
N1_cr_upperSide = 6.0;
eNformulation = 1;

% Parameters to detect BL edge:
eps0 = 0.01;
eps1 = 0.1;

numPointsPerElem_x = 1;
len_y = 0.06;                % Non-dimensionalized with respect to the chord

numExtrusionLayers = length(unique(mesh.p(:,3))) - 1;
extrusionLength = max(mesh.p(:,3)) - min(mesh.p(:,3));

if rem(timeStepEnd - timeStepStart, timeStepFreq) ~= 0; error('timeStepFreq not consistent with timeStepStart and timeStepEnd.'); end
if numPointsPerElem_x > 1; error('When we order the points in sField, xyCoord, etc. it may be wrong when numPointsPerElem_x > 1'); end
if mesh.nd ~= 3; error('LES post-processing only available for 3D simulations.'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% COMPUTATION STARTS HERE %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

torder = app.torder;
porder = mesh.porder;       % Polynomial order of the numerical approximation (in the master element)
pgauss = integralAccuracy*porder;          % Polynomial order up to which quadrature rule provides exact result
elemtype = mesh.elemtype;
plocvl = master.plocvl;
plocfc = master.plocfc;
npf = master.npf;
nfe = mesh.nfe;
nd = mesh.nd;
perm = mesh.perm;
dgnodes = mesh.dgnodes;
nc = app.nc;
npv = master.npv;
ne = mesh.ne;
nf = mesh.nf;
nch = app.nch;
pMesh = mesh.p;
tMesh = mesh.t;
fMesh = mesh.f;
t2fMesh = mesh.t2f;
gam = app.arg{1};
Re_inf = app.arg{3};
dt = app.dt(end);

dz = extrusionLength / numExtrusionLayers;
numTimeSteps = (timeStepEnd - timeStepStart) / timeStepFreq + 1;
if rem(numTimeSteps,2) == 1; warning('Not sure if FFT implementation will work fine with odd number of time-steps.'); end
samplingFreq = 1/dt;
frequencies = samplingFreq*(0:(numTimeSteps/2))/numTimeSteps;

for ii=0:porder
    pp_vl{ii+1} = jacobi(ii,0,0);
    pp_vl{ii+1} = pp_vl{ii+1}*sqrt(2*ii+1);
    dpp_vl{ii+1} = polyder(pp_vl{ii+1});
    pp_fc{ii+1} = jacobi(ii,0,0);
    pp_fc{ii+1} = pp_fc{ii+1}*sqrt(2*ii+1);
    dpp_fc{ii+1} = polyder(pp_fc{ii+1});
end

% Compute Alocvl and Alocfc
if elemtype == 0
    Alocvl = koornwinder(plocvl,porder,pp_vl,dpp_vl);
    Alocfc = koornwinder(plocfc,porder,pp_fc,dpp_fc);
else
    Alocvl = tensorproduct(plocvl,porder,pp_vl,dpp_vl);
    Alocfc = tensorproduct(plocfc,porder,pp_fc,dpp_fc);
end
Alocvl_inv = inv(Alocvl);
Alocfc_inv = inv(Alocfc);
clear Alocvl Alocfc

% Compute numPointsPerElem
if elemtype == 0
    aux_numPointsPerElem = tril(ones(numPointsPerElem_x,numPointsPerElem_x),0);
    numPointsPerElem = sum(aux_numPointsPerElem(:));
elseif elemtype == 1
    numPointsPerElem = numPointsPerElem_x*numPointsPerElem_x;
end
clear aux_numPointsPerElem


% Get faces, elements and relative face index next to the airfoil
facesOnAirfoil = find(mesh.f(:,end) == -ib);       % Faces on the airfoil surface
numFacesOnAirfoil = length(facesOnAirfoil);
numFacesOnAirfoilPerExtrusionLayer = numFacesOnAirfoil / numExtrusionLayers;
numElemPerExtrusionLayer = floor(ne / numExtrusionLayers + 0.01);
elemOnAirfoil = mesh.f(facesOnAirfoil,end-1);         % Elements next to the airfoil
numElemOnAirfoil = length(elemOnAirfoil);
numElemOnAirfoilPerExtrLayer = numElemOnAirfoil / numExtrusionLayers;

if rem(numFacesOnAirfoil,numExtrusionLayers) ~= 0; error('Number of faces on airfoil not compatible with number of extrusion layers.'); end
if rem(ne,numExtrusionLayers) ~= 0; error('Number of elements not compatible with number of extrusion layers.'); end

if min(elemOnAirfoil < 1) || max(elemOnAirfoil) > ne
    error('Something wrong.');
end

indexFaceOnAirfoil = zeros(numFacesOnAirfoil,1);
for i=1:numFacesOnAirfoil
    indexFaceOnAirfoil(i) = find(mesh.t2f(elemOnAirfoil(i),:) == facesOnAirfoil(i));     % Local index of the face on the airfoil
    
    if length(indexFaceOnAirfoil(i)) ~= 1
        error('Error No. 1');
    end
end


% Get faces on airfoil and first extrusion layer
verticesOnAirfoil = mesh.f(facesOnAirfoil,1:end-2);
zMin = min(mesh.p(verticesOnAirfoil(:),3));
[facesOnAirfoilFirstExtrLayer_aux, ~] = find(reshape(mesh.p(verticesOnAirfoil(:),3),size(verticesOnAirfoil)) < zMin+1.0e-8);
facesOnAirfoilFirstExtrLayer_aux = unique(facesOnAirfoilFirstExtrLayer_aux(:));
facesOnAirfoilFirstExtrLayer = facesOnAirfoil(facesOnAirfoilFirstExtrLayer_aux);
numFacesOnAirfoilFirstExtrLayer = length(facesOnAirfoilFirstExtrLayer);

if numFacesOnAirfoilFirstExtrLayer * numExtrusionLayers ~= numFacesOnAirfoil
    error('Number of faces on airfoil and first extrusion layer is not correct.');
end

% Map from facesOnAirfoil to facesOnAirfoilFirstExtrLayer
mapping_tmp = zeros(max(facesOnAirfoilFirstExtrLayer_aux),1);
mapping_tmp(facesOnAirfoilFirstExtrLayer_aux) = 1:numFacesOnAirfoilFirstExtrLayer;
mappingToFaceOn1stExtrLayer = zeros(numFacesOnAirfoil,1);
mappingToFaceOn1stExtrLayer_aux = zeros(numFacesOnAirfoil,1);
for i=1:numFacesOnAirfoil
    verticesOnNthLayer = mesh.p(mesh.f(facesOnAirfoil(i),1:end-2),:);
    for j=1:numFacesOnAirfoilFirstExtrLayer
        verticesOn1stLayer = mesh.p(mesh.f(facesOnAirfoilFirstExtrLayer(j),1:end-2),:);
        if all(sqrt(sum((sortrows(verticesOnNthLayer(:,1:2))-sortrows(verticesOn1stLayer(:,1:2))).^2,2)) < 1e-7) %max(sqrt(sum((verticesOnNthLayer(:,1:2)-verticesOn1stLayer(:,1:2)).^2,2)) < 1e-7)
            mappingToFaceOn1stExtrLayer(i) = facesOnAirfoilFirstExtrLayer(j);
            mappingToFaceOn1stExtrLayer_aux(i) = facesOnAirfoilFirstExtrLayer_aux(j);
            break;
        end
    end
end
mappingToFaceOn1stExtrLayer_aux_aux = mapping_tmp(mappingToFaceOn1stExtrLayer_aux);
if min(mappingToFaceOn1stExtrLayer) == 0; error('Face on airfoil and Nth extrusion layer was not match to face on first layer.'); end
if min(mappingToFaceOn1stExtrLayer_aux_aux) == 0; error('Face on airfoil and Nth extrusion layer was not match to face on first layer.'); end
clear mapping_tmp


% Mapping from elements in Nth layer to elements in 1st layer, and vice versa:
elemEquivToElemIn1stLayer = zeros(numExtrusionLayers,numElemPerExtrusionLayer);
elemToElemIn1stLayer = zeros(ne,1);
for i=1:numElemPerExtrusionLayer
    elemStart = i;
    elemEnd = i + (numExtrusionLayers-1)*numElemPerExtrusionLayer;
    elemEquivToElemIn1stLayer(:,i) = elemStart:numElemPerExtrusionLayer:elemEnd;
    elemToElemIn1stLayer(elemEquivToElemIn1stLayer(:,i)) = i;
end
if min(elemEquivToElemIn1stLayer(:)) < 1 || min(elemToElemIn1stLayer) < 1; error('Mapping from elements in Nth layer to elements in 1st layer went wrong.'); end
if max(elemEquivToElemIn1stLayer(:)) < ne || max(elemToElemIn1stLayer) < numElemPerExtrusionLayer; error('Mapping from elements in Nth layer to elements in 1st layer went wrong.'); end
clear elemStart elemEnd



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Courant number %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

masterFaceCentroid = getMasterCentroid(elemtype,nd-1);
masterFaceCentroid = reshape(masterFaceCentroid,[1,nd-1]);
shap_f = mkshape(porder, plocfc, masterFaceCentroid, elemtype, Alocfc_inv, pp_fc, dpp_fc);
if elemtype == 0
    hField = zeros(nfe,ne);          % Heights
    for element=1:ne
        for i=1:nfe
            vertex = mesh.p(mesh.t(element,i),:);
            vertex = vertex(:);
            midPointOppFace = squeeze(shap_f(:,1,1))' * squeeze(dgnodes(perm(:,i),1:3,element));
            midPointOppFace = midPointOppFace(:);
            hField(i,element) = norm(vertex-midPointOppFace);
        end
    end
elseif elemtype == 1
    hField = zeros(nfe/2,ne);        % Distance between parallel faces
    for element=1:ne
        i1 = -1;
        for i=1:nfe/2
            i1 = i1 + 2;
            i2 = i1 + 1;
            centroidFace1 = squeeze(shap_f(:,1,1))' * squeeze(dgnodes(perm(:,i1),1:3,element));
            centroidFace2 = squeeze(shap_f(:,1,1))' * squeeze(dgnodes(perm(:,i2),1:3,element));
            hField(i,element) = norm(centroidFace1-centroidFace2);
        end
    end
else
    error('elemtype not implemented.');
end
hMin = min(hField,[],1);     % Field of smallest height of each element
hMax = max(hField,[],1);     % Field of largest height of each element
hMinMin = min(hMin);    % For CFL condition
hMaxMax = max(hMax);

CourantNumber = dt / hMinMin;
disp('**********************************');
disp(['Courant number is ', num2str(CourantNumber)]);
disp(['Spanwise Courant number is ', num2str(dt/dz)]);
disp('**********************************');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute field of normals and coordinates of points on airfoil %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nField_tmp = cell(numFacesOnAirfoilFirstExtrLayer);
xFieldSurface = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,nd);
mesh_ = mesh;
for i=1:numFacesOnAirfoilFirstExtrLayer
    nField_tmp{i} = zeros(numPointsPerElem,nd);
    
    [coordPointsMaster2d, ~] = getCoordPointsMaster(mesh_, indexFaceOnAirfoil(facesOnAirfoilFirstExtrLayer_aux(i)), numPointsPerElem_x);
    shap2d = mkshape(porder, plocfc, coordPointsMaster2d, elemtype, Alocfc_inv, pp_fc, dpp_fc);
    dgnodesFace = dgnodes(perm(:,indexFaceOnAirfoil(facesOnAirfoilFirstExtrLayer_aux(i))),1:3,elemOnAirfoil(facesOnAirfoilFirstExtrLayer_aux(i)));
    
    % Compute normal to airfoil
    t1 = squeeze(shap2d(:,:,2))' * dgnodesFace;
    t2 = squeeze(shap2d(:,:,3))' * dgnodesFace;
    for l=1:numPointsPerElem
        n = - cross(t1(l,:),t2(l,:));       % Note: The minus sign ensures inwards-pointing normal
        nField_tmp{i}(l,:) = n/norm(n);
    end
    if max(abs(nField_tmp{i}(:,3))) > 1.0e-8; error('Normal vector with non-zero z-component.'); end
    
    % Compute coordinates of points on surface
    xFieldSurface_tmp = squeeze(shap2d(:,:,1))' * dgnodesFace;
    xFieldSurface(i,:,:) = reshape(xFieldSurface_tmp,[numPointsPerElem,nd]);
end

nField = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,nd);
for i=1:numFacesOnAirfoilFirstExtrLayer
    nField(i,:,:) = nField_tmp{i}(:,:);
end

clear mesh_



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get field of unique (x,y)-coordinates of points on airfoil %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snap = 1.0e-8;
xyCoord = round(xFieldSurface(:,:,1:2)/snap)*snap;
xyCoord = reshape(xyCoord,[],2);
[~,I] = unique(xyCoord,'rows');
I = sort(I);
xyCoordUnique = xyCoord(I,:);
[~,mapping1stLayerPoints2Uniq1stLayerPoints] = ismember(xyCoord, xyCoordUnique,'rows');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get field of s-coordinates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note : s direction connects LE and TE
sField = get_s_coord(LEcoord,alpha,reshape(xFieldSurface,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem,nd]));

sField_tmp = zeros(size(xyCoordUnique,1),1);
numPointsZ = zeros(size(xyCoordUnique,1),1);
if length(mapping1stLayerPoints2Uniq1stLayerPoints) ~= numFacesOnAirfoilFirstExtrLayer*numPointsPerElem_x
    warning('Number of sampling directions is unexpected.');
end
for i=1:length(mapping1stLayerPoints2Uniq1stLayerPoints)
    sField_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = sField_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + sField(i);
    numPointsZ(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = numPointsZ(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + 1;
end
if min(numPointsZ) < 1; error('Some points were not mapped properly.'); end
sField = sField_tmp ./ numPointsZ;

clear shap2d dgnodesFace t1 t2 coordPointsMaster2d nField_tmp mesh_ sField_tmp numPointsZ



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get indices in lower and upper side of airfoil %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% warning('Current code assumes some ordering of faces on airfoil. This may not work for non-c meshes.');
% [~,iLE] = min(sField(:));
% iLowerSide = 1:iLE;
% iUpperSide = iLE:length(sField);
% [sField(iLowerSide), iLowerSidePermIndices] = sort(sField(iLowerSide));
% sField(iLowerSide) = flip(sField(iLowerSide));
% iLowerSidePermIndices = flip(iLowerSidePermIndices);
% [sField(iUpperSide), iUpperSidePermIndices] = sort(sField(iUpperSide));
% sFieldLowerSide = sField(iLowerSide);
% sFieldUpperSide = sField(iUpperSide);
% iLowerSide = iLowerSide(iLowerSidePermIndices);
% iUpperSide = iUpperSide(iUpperSidePermIndices);

%%% T106C:
warning('ONLY VALID FOR T106C');
[~,iLE] = min(sField(:));
[~,iTE] = max(sField(:));
iLowerSide = [iLE:length(sField), 1:(iTE-1)];
iUpperSide = iTE:iLE;
[sField(iLowerSide), iLowerSidePermIndices] = sort(sField(iLowerSide));
sField(iLowerSide) = flip(sField(iLowerSide));
iLowerSidePermIndices = flip(iLowerSidePermIndices);
[sField(iUpperSide), iUpperSidePermIndices] = sort(sField(iUpperSide));
sFieldLowerSide = sField(iLowerSide);
sFieldUpperSide = sField(iUpperSide);
iLowerSide = iLowerSide(iLowerSidePermIndices);
iUpperSide = iUpperSide(iUpperSidePermIndices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get indices of points on airfoil associated to microphones %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, index relative to unique numbering of 2D points on surface
% numMicrophones = 6;
iUnique25PctLower = find(sFieldLowerSide > 0.25, 1, 'last');
iUnique50PctLower = find(sFieldLowerSide > 0.50, 1, 'last');
iUnique75PctLower = find(sFieldLowerSide > 0.75, 1, 'last');
iUnique25PctUpper = length(iLowerSide) + find(sFieldUpperSide > 0.25, 1, 'first');
iUnique50PctUpper = length(iLowerSide) + find(sFieldUpperSide > 0.50, 1, 'first');
iUnique75PctUpper = length(iLowerSide) + find(sFieldUpperSide > 0.75, 1, 'first');
[uniqueSfieldLower,iUniqueSfieldLower,~] = unique(sFieldLowerSide);
[uniqueSfieldUpper,iUniqueSfieldUpper,~] = unique(sFieldUpperSide);
iUniqueSfieldUpper = length(iLowerSide) - 1 + iUniqueSfieldUpper;

% % Then, index relative to non-unique numbering of first-layer points on surface
% xFieldSurface = reshape(xFieldSurface,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem,nd]);
% i25PctLower = find(mapping1stLayerPoints2Uniq1stLayerPoints == iUnique25PctLower);
% if isempty(i25PctLower); error('No candidate nodes detected.'); end
% [~,index] = min(xFieldSurface(i25PctLower,3));
% i25PctLower = i25PctLower(index);
% i50PctLower = find(mapping1stLayerPoints2Uniq1stLayerPoints == iUnique50PctLower);
% if isempty(i50PctLower); error('No candidate nodes detected.'); end
% [~,index] = min(xFieldSurface(i50PctLower,3));
% i50PctLower = i50PctLower(index);
% i75PctLower = find(mapping1stLayerPoints2Uniq1stLayerPoints == iUnique75PctLower);
% if isempty(i75PctLower); error('No candidate nodes detected.'); end
% [~,index] = min(xFieldSurface(i75PctLower,3));
% i75PctLower = i75PctLower(index);
% i25PctUpper = find(mapping1stLayerPoints2Uniq1stLayerPoints == iUnique25PctUpper);
% if isempty(i25PctUpper); error('No candidate nodes detected.'); end
% [~,index] = min(xFieldSurface(i25PctUpper,3));
% i25PctUpper = i25PctUpper(index);
% i50PctUpper = find(mapping1stLayerPoints2Uniq1stLayerPoints == iUnique50PctUpper);
% if isempty(i50PctUpper); error('No candidate nodes detected.'); end
% [~,index] = min(xFieldSurface(i50PctUpper,3));
% i50PctUpper = i50PctUpper(index);
% i75PctUpper = find(mapping1stLayerPoints2Uniq1stLayerPoints == iUnique75PctUpper);
% if isempty(i75PctUpper); error('No candidate nodes detected.'); end
% [~,index] = min(xFieldSurface(i75PctUpper,3));
% i75PctUpper = i75PctUpper(index);
% surfacePoints4Microphones = [i25PctLower, i50PctLower, i75PctLower, i25PctUpper, i50PctUpper, i75PctUpper];
% xFieldSurface = reshape(xFieldSurface,[numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,nd]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determine integration points for each element %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shapg = cell(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
UtoLineIntegralU = cell(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
elemInNormalDir = cell(numFacesOnAirfoil, numPointsPerElem);
lenInElem = cell(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
numElemInNormalDir = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);

[gaussPoints1d,gaussWeights1d] = gaussquad(pgauss,1,0);
numGaussPoints = length(gaussPoints1d);
gaussWeights1d = reshape(gaussWeights1d,[1,numGaussPoints]);        % Make row vector
gaussPoints1d = reshape(gaussPoints1d,[numGaussPoints,1]);          % Make column vector. 1D Gauss points in interval [0,1] . Size: [numGaussPoints, 1]
points1d = [0; gaussPoints1d; 1];
numPointsNormalDirPerElem = length(points1d);

% First, compute first extrusion layer:
for i=1:numFacesOnAirfoilFirstExtrLayer
    disp(['Face No. ',num2str(i),' / ', num2str(numFacesOnAirfoilPerExtrusionLayer)]);
    for j=1:numPointsPerElem
        numElements = 0;
        
        shapg{i,j} = {};
        UtoLineIntegralU{i,j} = {};
        elemInNormalDir{facesOnAirfoilFirstExtrLayer_aux(i),j} = {};
        lenInElem{i,j} = {};

        totalLength = 0;
        normal = nField(i,j,:);
        normal = normal(:);
        initialPoint = reshape(xFieldSurface(i,j,:),[1,nd]);
        currentLocalFace = indexFaceOnAirfoil(facesOnAirfoilFirstExtrLayer_aux(i));
        currentElem = elemOnAirfoil(facesOnAirfoilFirstExtrLayer_aux(i));
        
        while true
            numElements = numElements + 1;
            elemInNormalDir{facesOnAirfoilFirstExtrLayer_aux(i),j}{numElements} = currentElem;
            
            % Compute what face the normal with intersect with
            [nextFace_local, lenInElem{i,j}{numElements}, finalPoint] = getNextFace(normal, initialPoint, ...
                currentLocalFace, pMesh(tMesh(currentElem,:),:),porder,plocfc,dgnodes(:,1:3,currentElem), ...
                perm,elemtype,npf,Alocfc_inv, pp_fc, dpp_fc);
            initialPoint = reshape(initialPoint,[1,nd]);
            finalPoint = reshape(finalPoint,[1,nd]);
            nextFace = t2fMesh(currentElem,nextFace_local);
            
            % Compute integration points in the element
            UtoLineIntegralU{i,j}{numElements} = zeros(numPointsNormalDirPerElem, npv);
            pointsInElem = repmat(initialPoint,[numPointsNormalDirPerElem,1]) + points1d*(finalPoint-initialPoint);
            pointsInElem_master = mapBack(squeeze(dgnodes(:,1:3,currentElem)), pointsInElem, porder, elemtype, plocvl, npv, nd, Alocvl_inv, pp_vl, dpp_vl);
            shapg{i,j}{numElements} = mkshape(porder, plocvl, pointsInElem_master, elemtype, Alocvl_inv, pp_vl, dpp_vl);
            
            gaussNormalDir = kron(pointsInElem(1,:),ones(numGaussPoints*(numPointsNormalDirPerElem-1),1)) + kron(pointsInElem(2:end,:) - repmat(pointsInElem(1,:),[numPointsNormalDirPerElem-1,1]), gaussPoints1d);
            gaussNormalDirMaster = mapBack(squeeze(dgnodes(:,1:3,currentElem)), gaussNormalDir, porder, elemtype, plocvl, npv, nd, Alocvl_inv, pp_vl, dpp_vl);
            % Compute "shape" matrix that maps the solution at the nodes to
            % the solution at the desired Gauss points
            tmp22 = mkshape(porder, plocvl, gaussNormalDirMaster, elemtype, Alocvl_inv, pp_vl, dpp_vl);
            
            intLenDotW = lenInElem{i,j}{numElements} * points1d(2:end) * gaussWeights1d;
            for lll=2:numPointsNormalDirPerElem
                init1 = (lll-2)*numGaussPoints + 1;
                endd1 = init1 + numGaussPoints - 1;
                UtoLineIntegralU{i,j}{numElements}(lll, :) = intLenDotW(lll-1,:) * squeeze(tmp22(:,init1:endd1,1))';
            end

            % Determine next element
            nextElem = fMesh(nextFace,end-1:end);
            nextElem = setdiff(nextElem,currentElem);
            if length(nextElem) ~= 1; error('Error No. 2'); end
            if nextElem <= 0; error('Outside of problem domain'); end
            
            nextLocalFace = find(t2fMesh(nextElem,:) == nextFace);
            if length(nextLocalFace) ~= 1; error('More than one face detected for the next element.'); end
            
            totalLength = totalLength + lenInElem{i,j}{numElements};

            if totalLength > len_y; break; end

            currentElem = nextElem;
            currentLocalFace = nextLocalFace;
            initialPoint = finalPoint;
        end
        numElemInNormalDir(i,j) = numElements;
    end
end

% Then, extend results from first extrusion layer to all other extrusion layers:
for i=1:numFacesOnAirfoil
    iFirstLayer = mappingToFaceOn1stExtrLayer_aux_aux(i);
    
    verticesOnFace = mesh.f(facesOnAirfoil(i),1:end-2);
    zCoordOfVerticesOnFace = mesh.p(verticesOnFace,3);
    zMin = min(zCoordOfVerticesOnFace);
    zMax = max(zCoordOfVerticesOnFace);
    if abs(zMax - zMin - dz) > 1e-8; error('dz in face not consistent with length of elements in extrusion direction.'); end
    layerNo = floor(zMax / dz + 0.01);
    if layerNo < 1 || layerNo > numExtrusionLayers
        error('layerNo has invalid value.');
    end
    
    for j=1:numPointsPerElem
        for l=1:length(elemInNormalDir{iFirstLayer,j})
            % Note: This in based on the fact that the elements in the
            % first layer are ordered first, and there is an easy mapping
            % from elements in the Nth layer to elements in the 1st layer.
            elemInNormalDir{i,j}{l} = elemInNormalDir{iFirstLayer,j}{l} + floor((layerNo-1)*numElemPerExtrusionLayer + 0.01);
        end
    end
end
maxElemInNormalDir = max(numElemInNormalDir(:));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute time average UDG %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UDG_tAvg = zeros(npv, nc, ne);
UDG_Inst = zeros(npv, nc, ne);
% if strcmp(mesh.hybrid,'hdg')
%     UH = zeros(nch, npf, nf);
% elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
%     UH = zeros(nch, mesh.ndh);
% end

timeStepCounter = 0;
for m=timeStepStart:timeStepFreq:timeStepEnd
    timeStepCounter = timeStepCounter + 1;
    
    % Read solution file
    for i=1:numProc
        fileID = fopen([fileName, '_t', num2str(m), '_np', num2str(i-1), '.bin'], 'r');
        data = fread(fileID,'double');
        fclose(fileID);
        
        nElemIn = sum(dmd{i}.elempartpts(1:2));
%         nEntIn = sum(dmd{i}.entpartpts(1:2)); 
        indElem = dmd{i}.elempart(1:nElemIn);               
%         indEnt = dmd{i}.entpart(1:nEntIn);
        n1 = npv*nc*nElemIn;        
        
        if strcmp(mesh.hybrid,'hdg')
            UDG_Inst(:,:,indElem) = reshape(data(1:n1),[npv nc nElemIn]);
%             UH(:,:,indEnt) = reshape(data(n1+1:end),[nch npf nEntIn]);
        elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
            UDG_Inst(:,:,indElem) = reshape(data(1:n1),[npv nc nElemIn]);
%             UH(:,indEnt) = reshape(data(n1+1:end),[nch nEntIn]);
        end
    end
    UDG_tAvg = UDG_tAvg + UDG_Inst;
end
UDG_tAvg = UDG_tAvg / numTimeSteps;
clear data UDG_Inst nElemIn nEntIn indElem n1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute time-spanwise average UDG %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UDG_s2tAvg = zeros(npv, nc, numElemPerExtrusionLayer);
for i=1:numElemPerExtrusionLayer
    UDG_s2tAvg(:,:,i) = sum(UDG_tAvg(:,:,elemEquivToElemIn1stLayer(:,i)),3) / numExtrusionLayers;
end
r_s2tAvg = UDG_s2tAvg(:,1,:);
uv_s2tAvg = UDG_s2tAvg(:,2,:)./UDG_s2tAvg(:,1,:);
vv_s2tAvg = UDG_s2tAvg(:,3,:)./UDG_s2tAvg(:,1,:);
wv_s2tAvg = UDG_s2tAvg(:,4,:)./UDG_s2tAvg(:,1,:);
E_s2tAvg = UDG_s2tAvg(:,5,:)./UDG_s2tAvg(:,1,:);
p_s2tAvg = (gam-1)*r_s2tAvg.*(E_s2tAvg-0.5*(uv_s2tAvg.^2+vv_s2tAvg.^2+wv_s2tAvg.^2));
clear UDG_tAvg



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute [time-spanwise] (s2,t)-average boundary layer parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem);
deltaS_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem);
theta_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem);
thetaS_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem);
H_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem);
HS_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem);
u_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem, nd);
uMag_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
u1_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
u2_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
u3_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
uS_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem, nd);
uSmag_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
uS1_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
uS2_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
uS3_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
p_InNormal_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
ue_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, nd);
ueMag_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
ueS_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, nd);
ueSmag_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
pe_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
s1_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, nd);
s2_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem, nd);
ne_s2tAvg = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
h_n = zeros(numFacesOnAirfoilFirstExtrLayer);
u_tau = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
u1Plus = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
Du1Plus = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
yPlus = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
nDist_InNormal = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);
noElemInBL = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);

r_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
uv_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
ru_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
vv_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rv_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
wv_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rw_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
E_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rE_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
p_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
r_x_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
ru_x_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rv_x_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rw_x_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
r_y_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
ru_y_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rv_y_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rw_y_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
r_z_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
ru_z_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rv_z_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
rw_z_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
uv_x_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
vv_x_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
wv_x_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
uv_y_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
vv_y_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
wv_y_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
uv_z_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
vv_z_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
wv_z_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
vort_x_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,nd);
Dvort_xDn_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,nd);
vortCrossN_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,nd);
vortMag_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);
DvortDnMag_InElem_s2tAvg = zeros(maxElemInNormalDir*npv,1);

% microphone2nodeInElem = zeros(numMicrophones,1);
% microphone2elem = zeros(numMicrophones,1);

shap_plocvl = mkshape(porder, plocvl, plocvl, elemtype, Alocvl_inv, pp_vl, dpp_vl);
shap_plocvl = shap_plocvl(:,:,2:end);

for i=1:numFacesOnAirfoilFirstExtrLayer
    disp(['Face No. ', num2str(i), ' / ', num2str(numFacesOnAirfoilFirstExtrLayer)]);
    
    % Compute h_n:
    faceIndex = indexFaceOnAirfoil(facesOnAirfoilFirstExtrLayer_aux(i));
    elementNearFoil = elemToElemIn1stLayer(elemInNormalDir{facesOnAirfoilFirstExtrLayer_aux(i),1}{1});
    if elemtype == 0
        vertex = mesh.p(mesh.t(elementNearFoil,faceIndex),:);
        vertex = vertex(:);
        midPointOppFace = squeeze(shap_f(:,1,1))' * squeeze(dgnodes(perm(:,faceIndex),1:3,elementNearFoil));
        midPointOppFace = midPointOppFace(:);
        h_n(i) = abs(dot(vertex-midPointOppFace,squeeze(nField(i,j,:))));
    elseif elemtype == 1
        i1 = faceIndex;
        if rem(i1,2); i2 = i1 + 1; else i2 = i1 - 1; end
        centroidFace1 = squeeze(shap_f(:,1,1))' * squeeze(dgnodes(perm(:,i1),1:3,elementNearFoil));
        centroidFace2 = squeeze(shap_f(:,1,1))' * squeeze(dgnodes(perm(:,i2),1:3,elementNearFoil));
        h_n(i) = abs(dot(centroidFace1-centroidFace2,squeeze(nField(i,j,:))));
    end
    
    for j=1:numPointsPerElem
        currentNormal = nField(i,j,:);
        currentNormal = currentNormal(:);
        
        BLedgeFound = 0;
        
        % Compute BL velocities:
        for l=1:numElemInNormalDir(i,j)
            currentElem = elemInNormalDir{facesOnAirfoilFirstExtrLayer_aux(i),j}{l};
            currentElem1stLayer = elemToElemIn1stLayer(currentElem);
            
            init = (l-1) * npv + 1;
            endd = init + npv - 1;

            r_InElem_s2tAvg(init:endd) = UDG_s2tAvg(:,1,currentElem1stLayer);
            ru_InElem_s2tAvg(init:endd) = UDG_s2tAvg(:,2,currentElem1stLayer);
            rv_InElem_s2tAvg(init:endd) = UDG_s2tAvg(:,3,currentElem1stLayer);
            rw_InElem_s2tAvg(init:endd) = UDG_s2tAvg(:,4,currentElem1stLayer);
            rE_InElem_s2tAvg(init:endd) = UDG_s2tAvg(:,5,currentElem1stLayer);
            r_x_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,6,currentElem1stLayer);
            ru_x_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,7,currentElem1stLayer);
            rv_x_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,8,currentElem1stLayer);
            rw_x_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,9,currentElem1stLayer);
            r_y_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,11,currentElem1stLayer);
            ru_y_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,12,currentElem1stLayer);
            rv_y_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,13,currentElem1stLayer);
            rw_y_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,14,currentElem1stLayer);
            r_z_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,16,currentElem1stLayer);
            ru_z_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,17,currentElem1stLayer);
            rv_z_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,18,currentElem1stLayer);
            rw_z_InElem_s2tAvg(init:endd) = - UDG_s2tAvg(:,19,currentElem1stLayer);
            
            uv_InElem_s2tAvg(init:endd) = ru_InElem_s2tAvg(init:endd)./r_InElem_s2tAvg(init:endd);
            vv_InElem_s2tAvg(init:endd) = rv_InElem_s2tAvg(init:endd)./r_InElem_s2tAvg(init:endd);
            wv_InElem_s2tAvg(init:endd) = rw_InElem_s2tAvg(init:endd)./r_InElem_s2tAvg(init:endd);
            E_InElem_s2tAvg(init:endd) = rE_InElem_s2tAvg(init:endd)./r_InElem_s2tAvg(init:endd);
            p_InElem_s2tAvg(init:endd) = (gam-1)*r_InElem_s2tAvg(init:endd).*(E_InElem_s2tAvg(init:endd)-0.5*(uv_InElem_s2tAvg(init:endd).^2+vv_InElem_s2tAvg(init:endd).^2+wv_InElem_s2tAvg(init:endd).^2));
            
            uv_x_InElem_s2tAvg(init:endd) = (ru_x_InElem_s2tAvg(init:endd)-uv_InElem_s2tAvg(init:endd).*r_x_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            vv_x_InElem_s2tAvg(init:endd) = (rv_x_InElem_s2tAvg(init:endd)-vv_InElem_s2tAvg(init:endd).*r_x_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            wv_x_InElem_s2tAvg(init:endd) = (rw_x_InElem_s2tAvg(init:endd)-wv_InElem_s2tAvg(init:endd).*r_x_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            uv_y_InElem_s2tAvg(init:endd) = (ru_y_InElem_s2tAvg(init:endd)-uv_InElem_s2tAvg(init:endd).*r_y_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            vv_y_InElem_s2tAvg(init:endd) = (rv_y_InElem_s2tAvg(init:endd)-vv_InElem_s2tAvg(init:endd).*r_y_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            wv_y_InElem_s2tAvg(init:endd) = (rw_y_InElem_s2tAvg(init:endd)-wv_InElem_s2tAvg(init:endd).*r_y_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            uv_z_InElem_s2tAvg(init:endd) = (ru_z_InElem_s2tAvg(init:endd)-uv_InElem_s2tAvg(init:endd).*r_z_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            vv_z_InElem_s2tAvg(init:endd) = (rv_z_InElem_s2tAvg(init:endd)-vv_InElem_s2tAvg(init:endd).*r_z_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            wv_z_InElem_s2tAvg(init:endd) = (rw_z_InElem_s2tAvg(init:endd)-wv_InElem_s2tAvg(init:endd).*r_z_InElem_s2tAvg(init:endd))./r_InElem_s2tAvg(init:endd);
            
            initLength = 0;
            for ll=1:l-1
                initLength = initLength + lenInElem{i,j}{ll};
            end
            initLength = repmat(initLength,[numPointsNormalDirPerElem,1]);
            nDist_tmp = initLength + lenInElem{i,j}{l} * points1d;
            nDist_InNormal(i,j,l,:) = nDist_tmp;
            
            vort_x_InElem_s2tAvg(init:endd,:) = [wv_y_InElem_s2tAvg(init:endd)-vv_z_InElem_s2tAvg(init:endd), uv_z_InElem_s2tAvg(init:endd)-wv_x_InElem_s2tAvg(init:endd), vv_x_InElem_s2tAvg(init:endd)-uv_y_InElem_s2tAvg(init:endd)];
            
            term1 = reshape(shap_plocvl,[npv, npv*nd])' * squeeze(dgnodes(:,1:3,currentElem1stLayer));
            term1 = reshape(term1, [npv, nd, nd]);
            term1 = permute(term1, [2,3,1]);
            term2 = reshape(shap_plocvl,[npv, npv*nd])' * vort_x_InElem_s2tAvg(init:endd,1:3);
            term2 = reshape(term2, [npv, nd, nd]);
            term2 = permute(term2, [2,3,1]);
            for kk=1:npv
                Dvort_xDn_InElem_s2tAvg(init+kk-1,:) = currentNormal' * (squeeze(term1(:,:,kk))\squeeze(term2(:,:,kk)));
            end
            
            vortCrossN_InElem_s2tAvg(init:endd,:) = [vort_x_InElem_s2tAvg(init:endd,2)*currentNormal(3) - vort_x_InElem_s2tAvg(init:endd,3)*currentNormal(2), ...
                                       vort_x_InElem_s2tAvg(init:endd,3)*currentNormal(1) - vort_x_InElem_s2tAvg(init:endd,1)*currentNormal(3), ...
                                       vort_x_InElem_s2tAvg(init:endd,1)*currentNormal(2) - vort_x_InElem_s2tAvg(init:endd,2)*currentNormal(1)];
            
            vortMag_InElem_s2tAvg(init:endd) = sqrt(sum(vort_x_InElem_s2tAvg(init:endd,:).^2,2));
            DvortDnMag_InElem_s2tAvg(init:endd) = sqrt(sum(Dvort_xDn_InElem_s2tAvg(init:endd,:).^2,2));

            vortMag_aux = shapg{i,j}{l}(:,:,1)' * [vortMag_InElem_s2tAvg(init:endd), DvortDnMag_InElem_s2tAvg(init:endd)];
            vortMagIntPnts = vortMag_aux(:,1);
            DvortDnMagIntPnts = vortMag_aux(:,2);
            
            if l == 1
                uS_init = [0, 0, 0];
            else
                uS_init = uS_end;
            end

            aux_uS = UtoLineIntegralU{i,j}{l} * vortCrossN_InElem_s2tAvg(init:endd,:);
            uS_InNormal_s2tAvg(i,j,l,:,:) = repmat(uS_init,[numPointsNormalDirPerElem,1]) + aux_uS(:,:);
            uS_end = reshape(uS_InNormal_s2tAvg(i,j,l,end,:),[1,nd]);
            uSmag_InNormal_s2tAvg(i,j,l,:) = sqrt(sum(squeeze(uS_InNormal_s2tAvg(i,j,l,:,:)).^2,2));
            
            u_InNormal_s2tAvg(i,j,l,:,:) = shapg{i,j}{l}(:,:,1)' * [uv_InElem_s2tAvg(init:endd), vv_InElem_s2tAvg(init:endd), wv_InElem_s2tAvg(init:endd)];
            uMag_InNormal_s2tAvg(i,j,l,:) = sqrt(sum(squeeze(u_InNormal_s2tAvg(i,j,l,:,:)).^2,2));

            p_InNormal_s2tAvg(i,j,l,:) = shapg{i,j}{l}(:,:,1)' * p_InElem_s2tAvg(init:endd);
            
            criterion = (vortMagIntPnts .* nDist_tmp < eps0 * squeeze(uSmag_InNormal_s2tAvg(i,j,l,:))) .* (DvortDnMagIntPnts .* nDist_tmp.^2 < eps1 * squeeze(uSmag_InNormal_s2tAvg(i,j,l,:)));
%             criterion = (vortMagIntPnts .* nDist < eps0 * squeeze(uMag_stAvg(i,j,l,:))) .* (DvortDnMagIntPnts .* nDist.^2 < eps1 * squeeze(uMag_stAvg(i,j,l,:)));
            if all(criterion)
%                 [~,noPointsInLastElem] = max(criterion);
                noElemInBL(i,j) = l;
                BLedgeFound = 1;
                
                % Compute edge velocity and boundary layer thickness:
                ueS_s2tAvg(i,j,:) = uS_InNormal_s2tAvg(i,j,noElemInBL(i,j),numPointsNormalDirPerElem,:);
                ueSmag_s2tAvg(i,j) = norm(squeeze(ueS_s2tAvg(i,j,:)));
                ue_s2tAvg(i,j,:) = u_InNormal_s2tAvg(i,j,noElemInBL(i,j),numPointsNormalDirPerElem,:);
                ueMag_s2tAvg(i,j) = norm(squeeze(ue_s2tAvg(i,j,:)));
                ne_s2tAvg(i,j) = nDist_tmp(end);
                pe_s2tAvg(i,j) = p_InNormal_s2tAvg(i,j,noElemInBL(i,j),numPointsNormalDirPerElem);
                
                break;
            end
        end
        
        if BLedgeFound == 0;
            warning('Boundary layer edge was not found.');
            noElemInBL(i,j) = numElemInNormalDir(i,j);

            % Compute edge velocity (ue), edge pressure (pe) and boundary layer
            % thickness (ne)
            ueS_s2tAvg(i,j,:) = uS_InNormal_s2tAvg(i,j,noElemInBL(i,j),numPointsNormalDirPerElem,:);
            ueSmag_s2tAvg(i,j) = norm(squeeze(ueS_s2tAvg(i,j,:)));
            ue_s2tAvg(i,j,:) = u_InNormal_s2tAvg(i,j,noElemInBL(i,j),numPointsNormalDirPerElem,:);
            ueMag_s2tAvg(i,j) = norm(squeeze(ue_s2tAvg(i,j,:)));
            pe_s2tAvg(i,j) = p_InNormal_s2tAvg(i,j,noElemInBL(i,j),numPointsNormalDirPerElem);
            ne_s2tAvg(i,j) = nDist_tmp(end);
        end
        
%         % Assign element and node to microphone (if necessary)
%         index = find((j-1)*numFacesOnAirfoilFirstExtrLayer+i == surfacePoints4Microphones);
%         if length(index) > 1; error('Detection of location of BL in which microphones are located went wrong.'); end
%         if index
%             microphone2nodeInElem(index) = 1;     % TODO: Use a better approach
%             microphone2elem(index) = elemInNormalDir{facesOnAirfoilFirstExtrLayer_aux(i),j}{ceil(noElemInBL(i,j)/2)};
%             if microphone2elem(index) > numElemPerExtrusionLayer; error('Microphones are not located in first extrusion layer.'); end
%         end
        
        % Compute s1 and s2 unit vectors
        s1_s2tAvg(i,j,:) = ue_s2tAvg(i,j,:) / ueMag_s2tAvg(i,j);         % Components in x-coordinates
        s2_tmp = cross(squeeze(s1_s2tAvg(i,j,:)), currentNormal);
        s2_s2tAvg(i,j,:) = s2_tmp / norm(s2_tmp);          % Components in x-coordinates
        
        % Compute streamwise, normal and crossflow velocity components
        for l=1:noElemInBL(i,j)
            aux_uS = squeeze(uS_InNormal_s2tAvg(i,j,l,:,:)) * [squeeze(s1_s2tAvg(i,j,:)), squeeze(nField(i,j,:)), squeeze(s2_s2tAvg(i,j,:))];
            uS1_InNormal_s2tAvg(i,j,l,:) = aux_uS(:,1);
            uS2_InNormal_s2tAvg(i,j,l,:) = aux_uS(:,2);
            uS3_InNormal_s2tAvg(i,j,l,:) = aux_uS(:,3);
            aux_u = squeeze(u_InNormal_s2tAvg(i,j,l,:,:)) * [squeeze(s1_s2tAvg(i,j,:)), squeeze(nField(i,j,:)), squeeze(s2_s2tAvg(i,j,:))];
            u1_InNormal_s2tAvg(i,j,l,:) = aux_u(:,1);
            u2_InNormal_s2tAvg(i,j,l,:) = aux_u(:,2);
            u3_InNormal_s2tAvg(i,j,l,:) = aux_u(:,3);
        end
        
        % Compute y+:
        uv_x_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), uv_x_InElem_s2tAvg(1:npv));
        vv_x_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), vv_x_InElem_s2tAvg(1:npv));
        wv_x_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), wv_x_InElem_s2tAvg(1:npv));
        uv_y_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), uv_y_InElem_s2tAvg(1:npv));
        vv_y_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), vv_y_InElem_s2tAvg(1:npv));
        wv_y_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), wv_y_InElem_s2tAvg(1:npv));
        uv_z_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), uv_z_InElem_s2tAvg(1:npv));
        vv_z_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), vv_z_InElem_s2tAvg(1:npv));
        wv_z_wall_tmp = dot(squeeze(shapg{i,j}{1}(:,1,1)), wv_z_InElem_s2tAvg(1:npv));
        Dus1Dn_tmp = reshape(s1_s2tAvg(i,j,:),1,nd) * [uv_x_wall_tmp, uv_y_wall_tmp, uv_z_wall_tmp; vv_x_wall_tmp, vv_y_wall_tmp, vv_z_wall_tmp; wv_x_wall_tmp, wv_y_wall_tmp, wv_z_wall_tmp] * currentNormal(:);           % Here we assume DnDs1 is small
        u_tau(i,j) = sqrt((1 / Re_inf) * abs(Dus1Dn_tmp));     % Here we assume kinematic viscosity is constant
        u1Plus(i,j,:,:) = u1_InNormal_s2tAvg(i,j,:,:) / u_tau(i,j);
        Du1Plus(i,j,:,:) = (u1_InNormal_s2tAvg(i,j,:,:) - ue_s2tAvg(i,j,1)) / u_tau(i,j);
        yPlus(i,j) = Re_inf * h_n(i) * u_tau(i,j) / porder;     % Here we assume again kinematic viscosity is constant
        
        % Compute streamwise displacement and momentum thickness, and shape
        % parameter. What about crossflow thicknesses?
        for l=1:noElemInBL(i,j)
            aux_left = gaussWeights1d * lenInElem{i,j}{l}; % ones(1,numPointsNormalDirPerElem-1) * lenInElem{i,j}{l}/(numPointsNormalDirPerElem-1);
            
            % Displacement thickness
            aux_right = (1-squeeze(uS1_InNormal_s2tAvg(i,j,l,2:end-1))/ueSmag_s2tAvg(i,j));
            deltaS_s2tAvg(i,j) = deltaS_s2tAvg(i,j) + dot(aux_left,aux_right);
            aux_right = (1-squeeze(u1_InNormal_s2tAvg(i,j,l,2:end-1))/ueMag_s2tAvg(i,j));
            delta_s2tAvg(i,j) = delta_s2tAvg(i,j) + dot(aux_left,aux_right);
            
            % Momentum thickness
            aux_right = (1-squeeze(uS1_InNormal_s2tAvg(i,j,l,2:end-1))/ueSmag_s2tAvg(i,j)).*(squeeze(uS1_InNormal_s2tAvg(i,j,l,2:end-1))/ueSmag_s2tAvg(i,j));
            thetaS_s2tAvg(i,j) = thetaS_s2tAvg(i,j) + dot(aux_left,aux_right);
            aux_right = (1-squeeze(u1_InNormal_s2tAvg(i,j,l,2:end-1))/ueMag_s2tAvg(i,j)).*(squeeze(u1_InNormal_s2tAvg(i,j,l,2:end-1))/ueMag_s2tAvg(i,j));
            theta_s2tAvg(i,j) = theta_s2tAvg(i,j) + dot(aux_left,aux_right);
        end
        
        % Shape parameter
        HS_s2tAvg(i,j) = deltaS_s2tAvg(i,j) / thetaS_s2tAvg(i,j);
        H_s2tAvg(i,j) = delta_s2tAvg(i,j) / theta_s2tAvg(i,j);
    end
end
clear aux_uS aux_u

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute perturbation boundary layer parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UDG_Inst = zeros(npv, nc, ne);

u__InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem, nd);
uMag_InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
u1_InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
u2_InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
u3_InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
u_pert_InNormal = zeros(numFacesOnAirfoil,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem,nd,numTimeSteps);
u_pert_InNormal_s2tRMS = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem,nd);
uS_InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem, nd);
uSmag_InNormal = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
uS1_InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
uS2_InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
uS3_InNormal_Inst = zeros(numFacesOnAirfoil, numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem);
uS_pert_InNormal = zeros(numFacesOnAirfoil,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem,nd,numTimeSteps);
uS_pert_InNormal_s2tRMS = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem,nd);
pPert_InNormal = zeros(numFacesOnAirfoil,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem,numTimeSteps);
pPert_InNormal_s2tRMS = zeros(numFacesOnAirfoilFirstExtrLayer,numPointsPerElem,maxElemInNormalDir,numPointsNormalDirPerElem);

r_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
ru_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
rv_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
rw_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
rE_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
r_x_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
rv_x_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
rw_x_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
r_y_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
ru_y_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
rw_y_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
r_z_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
ru_z_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
rv_z_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
uv_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
vv_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
wv_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
E_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
p_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
vv_x_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
wv_x_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
uv_y_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
wv_y_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
uv_z_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
vv_z_InElem_Inst = zeros(maxElemInNormalDir*npv,1);
vort_x_InElem_Inst = zeros(maxElemInNormalDir*npv,3);
vortCrossN_InElem_Inst = zeros(maxElemInNormalDir*npv,3);
k_tAvg = zeros(npv,ne);

% pHistOnMicrophones = zeros(numTimeSteps,numMicrophones);
% kHistOnMicrophones = zeros(numTimeSteps,numMicrophones);

timeStepCounter = 0;
for m=timeStepStart:timeStepFreq:timeStepEnd
    timeStepCounter = timeStepCounter + 1;
    
    disp(['TIME-STEP NO. ', num2str(timeStepCounter), ' / ', num2str(numTimeSteps)]);
    
    % Read solution file
    for i=1:numProc
        fileID = fopen([fileName, '_t', num2str(m), '_np', num2str(i-1), '.bin'], 'r');
        data = fread(fileID,'double');
        fclose(fileID);
        
        nElemIn = sum(dmd{i}.elempartpts(1:2));
        indElem = dmd{i}.elempart(1:nElemIn);               
        n1 = npv*nc*nElemIn;        
        
        if strcmp(mesh.hybrid,'hdg')
            UDG_Inst(:,:,indElem) = reshape(data(1:n1),[npv nc nElemIn]);
        elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
            UDG_Inst(:,:,indElem) = reshape(data(1:n1),[npv nc nElemIn]);
        end
    end
    
    % Compute contribution to time-average turbulent kinetic energy:
    elements_tmp = repmat(1:numElemPerExtrusionLayer,1,numExtrusionLayers);
    k_tAvg(:,:) = k_tAvg(:,:) + squeeze(0.5*((UDG_Inst(:,2,:)./UDG_Inst(:,1,:) - UDG_s2tAvg(:,2,elements_tmp)./UDG_s2tAvg(:,1,elements_tmp)).^2 + ...
                                     (UDG_Inst(:,3,:)./UDG_Inst(:,1,:) - UDG_s2tAvg(:,3,elements_tmp)./UDG_s2tAvg(:,1,elements_tmp)).^2 + ...
                                     (UDG_Inst(:,4,:)./UDG_Inst(:,1,:) - UDG_s2tAvg(:,4,elements_tmp)./UDG_s2tAvg(:,1,elements_tmp)).^2));
    
%     % Get pressure and turbulent kinetic energy on microphones:
%     for i=1:numMicrophones
%         nodeInElem_tmp = microphone2nodeInElem(i);
%         elem_tmp = microphone2elem(i);
%         
%         r_tmp    = UDG_Inst(nodeInElem_tmp,1,elem_tmp);
%         uv_tmp   = UDG_Inst(nodeInElem_tmp,2,elem_tmp)./UDG_Inst(nodeInElem_tmp,1,elem_tmp);
%         vv_tmp   = UDG_Inst(nodeInElem_tmp,3,elem_tmp)./UDG_Inst(nodeInElem_tmp,1,elem_tmp);
%         wv_tmp   = UDG_Inst(nodeInElem_tmp,4,elem_tmp)./UDG_Inst(nodeInElem_tmp,1,elem_tmp);
%         E_tmp    = UDG_Inst(nodeInElem_tmp,5,elem_tmp)./UDG_Inst(nodeInElem_tmp,1,elem_tmp);
%         pCurrent_tmp = (gam-1)*r_tmp.*(E_tmp-0.5*(uv_tmp.^2+vv_tmp.^2+wv_tmp.^2));
%         
%         uvPert_tmp = uv_tmp - squeeze(uv_s2tAvg(nodeInElem_tmp,elem_tmp));
%         vvPert_tmp = vv_tmp - squeeze(vv_s2tAvg(nodeInElem_tmp,elem_tmp));
%         wvPert_tmp = wv_tmp - squeeze(wv_s2tAvg(nodeInElem_tmp,elem_tmp));
%         
%         kCurrent_tmp = 0.5 * (uvPert_tmp.^2 + vvPert_tmp.^2 + wvPert_tmp.^2);
%         
%         pHistOnMicrophones(timeStepCounter,i) = pHistOnMicrophones(timeStepCounter,i) + pCurrent_tmp;
%         kHistOnMicrophones(timeStepCounter,i) = kHistOnMicrophones(timeStepCounter,i) + kCurrent_tmp;
%     end
    
    for i=1:numFacesOnAirfoil
        if rem(i,200) == 1; disp(['Face No. ', num2str(i), ' / ', num2str(numFacesOnAirfoil)]); end
        
        for j=1:numPointsPerElem
            iFirstLayer = mappingToFaceOn1stExtrLayer_aux_aux(i);
            
            currentNormal = nField(iFirstLayer,j,:);
            currentNormal = currentNormal(:);

            for l=1:noElemInBL(iFirstLayer,j)
                currentElem = elemInNormalDir{i,j}{l};
                init = (l-1) * npv + 1;
                endd = init + npv - 1;

                r_InElem_Inst(init:endd) = UDG_Inst(:,1,currentElem);
                ru_InElem_Inst(init:endd) = UDG_Inst(:,2,currentElem);
                rv_InElem_Inst(init:endd) = UDG_Inst(:,3,currentElem);
                rw_InElem_Inst(init:endd) = UDG_Inst(:,4,currentElem);
                rE_InElem_Inst(init:endd) = UDG_Inst(:,5,currentElem);
                r_x_InElem_Inst(init:endd) = - UDG_Inst(:,6,currentElem);
                rv_x_InElem_Inst(init:endd) = - UDG_Inst(:,8,currentElem);
                rw_x_InElem_Inst(init:endd) = - UDG_Inst(:,9,currentElem);
                r_y_InElem_Inst(init:endd) = - UDG_Inst(:,11,currentElem);
                ru_y_InElem_Inst(init:endd) = - UDG_Inst(:,12,currentElem);
                rw_y_InElem_Inst(init:endd) = - UDG_Inst(:,14,currentElem);
                r_z_InElem_Inst(init:endd) = - UDG_Inst(:,16,currentElem);
                ru_z_InElem_Inst(init:endd) = - UDG_Inst(:,17,currentElem);
                rv_z_InElem_Inst(init:endd) = - UDG_Inst(:,18,currentElem);

                uv_InElem_Inst(init:endd) = ru_InElem_Inst(init:endd)./r_InElem_Inst(init:endd);
                vv_InElem_Inst(init:endd) = rv_InElem_Inst(init:endd)./r_InElem_Inst(init:endd);
                wv_InElem_Inst(init:endd) = rw_InElem_Inst(init:endd)./r_InElem_Inst(init:endd);
                E_InElem_Inst(init:endd) = rE_InElem_Inst(init:endd)./r_InElem_Inst(init:endd);
                p_InElem_Inst(init:endd) = (gam-1)*r_InElem_Inst(init:endd).*(E_InElem_Inst(init:endd)-0.5*(uv_InElem_Inst(init:endd).^2+vv_InElem_Inst(init:endd).^2+wv_InElem_Inst(init:endd).^2));
                vv_x_InElem_Inst(init:endd) = (rv_x_InElem_Inst(init:endd)-vv_InElem_Inst(init:endd).*r_x_InElem_Inst(init:endd))./r_InElem_Inst(init:endd);
                wv_x_InElem_Inst(init:endd) = (rw_x_InElem_Inst(init:endd)-wv_InElem_Inst(init:endd).*r_x_InElem_Inst(init:endd))./r_InElem_Inst(init:endd);
                uv_y_InElem_Inst(init:endd) = (ru_y_InElem_Inst(init:endd)-uv_InElem_Inst(init:endd).*r_y_InElem_Inst(init:endd))./r_InElem_Inst(init:endd);
                wv_y_InElem_Inst(init:endd) = (rw_y_InElem_Inst(init:endd)-wv_InElem_Inst(init:endd).*r_y_InElem_Inst(init:endd))./r_InElem_Inst(init:endd);
                uv_z_InElem_Inst(init:endd) = (ru_z_InElem_Inst(init:endd)-uv_InElem_Inst(init:endd).*r_z_InElem_Inst(init:endd))./r_InElem_Inst(init:endd);
                vv_z_InElem_Inst(init:endd) = (rv_z_InElem_Inst(init:endd)-vv_InElem_Inst(init:endd).*r_z_InElem_Inst(init:endd))./r_InElem_Inst(init:endd);
                
                vort_x_InElem_Inst(init:endd,:) = [wv_y_InElem_Inst(init:endd)-vv_z_InElem_Inst(init:endd), uv_z_InElem_Inst(init:endd)-wv_x_InElem_Inst(init:endd), vv_x_InElem_Inst(init:endd)-uv_y_InElem_Inst(init:endd)];

                vortCrossN_InElem_Inst(init:endd,:) = [vort_x_InElem_Inst(init:endd,2)*currentNormal(3) - vort_x_InElem_Inst(init:endd,3)*currentNormal(2), ...
                                           vort_x_InElem_Inst(init:endd,3)*currentNormal(1) - vort_x_InElem_Inst(init:endd,1)*currentNormal(3), ...
                                           vort_x_InElem_Inst(init:endd,1)*currentNormal(2) - vort_x_InElem_Inst(init:endd,2)*currentNormal(1)];

                if l == 1
                    uS_init = [0, 0, 0];
                else
                    uS_init = uS_end;
                end
                
                aux_uS = UtoLineIntegralU{iFirstLayer,j}{l} * vortCrossN_InElem_Inst(init:endd,:);
                uS_InNormal_Inst(i,j,l,:,:) = repmat(uS_init,[numPointsNormalDirPerElem,1]) + aux_uS(:,:);
                uS_end = reshape(uS_InNormal_Inst(i,j,l,end,:),[1,nd]);
                uSmag_InNormal(i,j,l,:) = sqrt(sum(squeeze(uS_InNormal_Inst(i,j,l,:,:)).^2,2));
                
                u__InNormal_Inst(i,j,l,:,:) = shapg{iFirstLayer,j}{l}(:,:,1)' * [uv_InElem_Inst(init:endd), vv_InElem_Inst(init:endd), wv_InElem_Inst(init:endd)];
                uMag_InNormal_Inst(i,j,l,:) = sqrt(sum(squeeze(u__InNormal_Inst(i,j,l,:,:)).^2,2));
                
                currentElem1stExtrLayer = rem(currentElem,numElemPerExtrusionLayer);
                if currentElem1stExtrLayer == 0; currentElem1stExtrLayer = numElemPerExtrusionLayer; end
                p_InNormal_s2tAvg = shapg{iFirstLayer,j}{l}(:,:,1)' * squeeze(p_s2tAvg(:,currentElem1stExtrLayer));
                p_InNormal_Inst = shapg{iFirstLayer,j}{l}(:,:,1)' * p_InElem_Inst(init:endd);
                pPert_InNormal(i,j,l,:,timeStepCounter) = p_InNormal_Inst - p_InNormal_s2tAvg;
            end

            % Compute instantaneous velocity and perturbation velocity:
            for l=1:noElemInBL(iFirstLayer,j)
                aux = squeeze(uS_InNormal_Inst(i,j,l,:,:)) * [squeeze(s1_s2tAvg(iFirstLayer,j,:)), squeeze(nField(iFirstLayer,j,:)), squeeze(s2_s2tAvg(iFirstLayer,j,:))];
                uS1_InNormal_Inst(i,j,l,:) = aux(:,1);
                uS2_InNormal_Inst(i,j,l,:) = aux(:,2);
                uS3_InNormal_Inst(i,j,l,:) = aux(:,3);
                
                uS_pert_InNormal(i,j,l,:,:,timeStepCounter) = [squeeze(uS1_InNormal_Inst(i,j,l,:)) - squeeze(uS1_InNormal_s2tAvg(iFirstLayer,j,l,:)), ...
                                                               squeeze(uS2_InNormal_Inst(i,j,l,:)) - squeeze(uS2_InNormal_s2tAvg(iFirstLayer,j,l,:)), ...
                                                               squeeze(uS3_InNormal_Inst(i,j,l,:)) - squeeze(uS3_InNormal_s2tAvg(iFirstLayer,j,l,:))];
                
                uS_pert_InNormal_s2tRMS(iFirstLayer,j,l,:,:) = uS_pert_InNormal_s2tRMS(iFirstLayer,j,l,:,:) + uS_pert_InNormal(i,j,l,:,:,timeStepCounter).^2;
                
                aux = squeeze(u__InNormal_Inst(i,j,l,:,:)) * [squeeze(s1_s2tAvg(iFirstLayer,j,:)), squeeze(nField(iFirstLayer,j,:)), squeeze(s2_s2tAvg(iFirstLayer,j,:))];
                u1_InNormal_Inst(i,j,l,:) = aux(:,1);
                u2_InNormal_Inst(i,j,l,:) = aux(:,2);
                u3_InNormal_Inst(i,j,l,:) = aux(:,3);
                
                u_pert_InNormal(i,j,l,:,:,timeStepCounter) = [squeeze(u1_InNormal_Inst(i,j,l,:)) - squeeze(u1_InNormal_s2tAvg(iFirstLayer,j,l,:)), ...
                                                              squeeze(u2_InNormal_Inst(i,j,l,:)) - squeeze(u2_InNormal_s2tAvg(iFirstLayer,j,l,:)), ...
                                                              squeeze(u3_InNormal_Inst(i,j,l,:)) - squeeze(u3_InNormal_s2tAvg(iFirstLayer,j,l,:))];
                
                u_pert_InNormal_s2tRMS(iFirstLayer,j,l,:,:) = u_pert_InNormal_s2tRMS(iFirstLayer,j,l,:,:) + u_pert_InNormal(i,j,l,:,:,timeStepCounter).^2;
                
                pPert_InNormal_s2tRMS(iFirstLayer,j,l,:) = pPert_InNormal_s2tRMS(iFirstLayer,j,l,:) + pPert_InNormal(i,j,l,:,timeStepCounter).^2;
            end
        end
    end
end
uS_pert_InNormal_s2tRMS = sqrt(uS_pert_InNormal_s2tRMS / (numTimeSteps*numExtrusionLayers));
u_pert_InNormal_s2tRMS = sqrt(u_pert_InNormal_s2tRMS / (numTimeSteps*numExtrusionLayers));
pPert_InNormal_s2tRMS = sqrt(pPert_InNormal_s2tRMS / (numTimeSteps*numExtrusionLayers));
% k_tAvg = k_tAvg / numTimeSteps;

% FFTpAtMicrophones = fft(pHistOnMicrophones);
% FFTpAtMicrophones = abs(FFTpAtMicrophones/numTimeSteps);
% FFTpAtMicrophones = FFTpAtMicrophones(1:numTimeSteps/2+1,:);
% FFTpAtMicrophones(2:end-1,:) = 2*FFTpAtMicrophones(2:end-1,:);
% 
% FFTkAtMicrophones = fft(kHistOnMicrophones);
% FFTkAtMicrophones = abs(FFTkAtMicrophones/numTimeSteps);
% FFTkAtMicrophones = FFTkAtMicrophones(1:numTimeSteps/2+1,:);
% FFTkAtMicrophones(2:end-1,:) = 2*FFTkAtMicrophones(2:end-1,:);

clear UDG_Inst u_ uMag u1 u2 u3 u_pert uS uSmag uS1 uS2 uS3 uS_pert



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute amplification of streamwise and crossflow perturbations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1 = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
An = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
A2 = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
AS1 = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
ASn = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
AS2 = zeros(numFacesOnAirfoilFirstExtrLayer, numPointsPerElem);
for i=1:numFacesOnAirfoilFirstExtrLayer
    for j=1:numPointsPerElem
        A_aux = [0, 0, 0];
        AS_aux = [0, 0, 0];
        
        for l=1:noElemInBL(i,j)
            aux_left = gaussWeights1d * lenInElem{i,j}{l};
            
            aux_right = squeeze(u_pert_InNormal_s2tRMS(i,j,l,2:end-1,:).^2);
            A_aux = A_aux + aux_left * aux_right;
            
            aux_right = squeeze(uS_pert_InNormal_s2tRMS(i,j,l,2:end-1,:).^2);
            AS_aux = AS_aux + aux_left * aux_right;
        end
        
        A_aux = sqrt(A_aux) ./ (ueMag_s2tAvg(i,j) * sqrt(ne_s2tAvg(i,j)));
        A1(i,j) = A_aux(1);
        An(i,j) = A_aux(2);
        A2(i,j) = A_aux(3);
        
        AS_aux = sqrt(AS_aux) ./ (ueSmag_s2tAvg(i,j) * sqrt(ne_s2tAvg(i,j)));
        AS1(i,j) = AS_aux(1);
        ASn(i,j) = AS_aux(2);
        AS2(i,j) = AS_aux(3);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Merge all points associated to the same (x,y) coordinates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_stAvg_tmp = zeros(size(xyCoordUnique,1),1);
theta_stAvg_tmp = zeros(size(xyCoordUnique,1),1);
H_stAvg_tmp = zeros(size(xyCoordUnique,1),1);
A1_tmp = zeros(size(xyCoordUnique,1),1);
An_tmp = zeros(size(xyCoordUnique,1),1);
A2_tmp = zeros(size(xyCoordUnique,1),1);
deltaS_stAvg_tmp = zeros(size(xyCoordUnique,1),1);
thetaS_stAvg_tmp = zeros(size(xyCoordUnique,1),1);
HS_stAvg_tmp = zeros(size(xyCoordUnique,1),1);
AS1_tmp = zeros(size(xyCoordUnique,1),1);
ASn_tmp = zeros(size(xyCoordUnique,1),1);
AS2_tmp = zeros(size(xyCoordUnique,1),1);
u_pert_s2tRMS_tmp = zeros(size(xyCoordUnique,1), maxElemInNormalDir, numPointsNormalDirPerElem, nd);
uS_pert_s2tRMS_tmp = zeros(size(xyCoordUnique,1), maxElemInNormalDir, numPointsNormalDirPerElem, nd);
u1_stAvg_tmp = zeros(size(xyCoordUnique,1), maxElemInNormalDir, numPointsNormalDirPerElem);
uS1_stAvg_tmp = zeros(size(xyCoordUnique,1), maxElemInNormalDir, numPointsNormalDirPerElem);
u1Plus_tmp = zeros(size(xyCoordUnique,1), maxElemInNormalDir, numPointsNormalDirPerElem);
Du1Plus_tmp = zeros(size(xyCoordUnique,1), maxElemInNormalDir, numPointsNormalDirPerElem);
pPert_s2tRMS_tmp = zeros(size(xyCoordUnique,1), maxElemInNormalDir, numPointsNormalDirPerElem);
nDist_InNormal_tmp = zeros(size(xyCoordUnique,1), maxElemInNormalDir, numPointsNormalDirPerElem);
pe_stAvg_tmp = zeros(size(xyCoordUnique,1),1);
ueSmag_s2tAvg_tmp = zeros(size(xyCoordUnique,1),1);
ueMag_s2tAvg_tmp = zeros(size(xyCoordUnique,1),1);
u_tau_tmp = zeros(size(xyCoordUnique,1),1);
yPlus_tmp = zeros(size(xyCoordUnique,1),1);
h_n_tmp = zeros(size(xyCoordUnique,1),1);
numPointsZ = zeros(size(xyCoordUnique,1),1);

u_pert_InNormal_s2tRMS = reshape(u_pert_InNormal_s2tRMS,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem, maxElemInNormalDir,numPointsNormalDirPerElem, nd]);
uS_pert_InNormal_s2tRMS = reshape(uS_pert_InNormal_s2tRMS,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem, maxElemInNormalDir,numPointsNormalDirPerElem, nd]);
u1_InNormal_s2tAvg = reshape(u1_InNormal_s2tAvg,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem]);
u1Plus = reshape(u1Plus,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem]);
Du1Plus = reshape(Du1Plus,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem]);
uS1_InNormal_s2tAvg = reshape(uS1_InNormal_s2tAvg,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem, maxElemInNormalDir, numPointsNormalDirPerElem]);
pPert_InNormal_s2tRMS = reshape(pPert_InNormal_s2tRMS,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem, maxElemInNormalDir,numPointsNormalDirPerElem]);
nDist_InNormal = reshape(nDist_InNormal,[numFacesOnAirfoilFirstExtrLayer*numPointsPerElem, maxElemInNormalDir,numPointsNormalDirPerElem]);
h_n = repmat(h_n(:), [numPointsPerElem,1]);
for i=1:length(mapping1stLayerPoints2Uniq1stLayerPoints)
    delta_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = delta_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + delta_s2tAvg(i);
    theta_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = theta_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + theta_s2tAvg(i);
    H_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = H_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + H_s2tAvg(i);
    A1_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = A1_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + A1(i);
    An_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = An_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + An(i);
    A2_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = A2_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + A2(i);
    deltaS_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = deltaS_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + deltaS_s2tAvg(i);
    thetaS_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = thetaS_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + thetaS_s2tAvg(i);
    HS_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = HS_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + HS_s2tAvg(i);
    AS1_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = AS1_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + AS1(i);
    ASn_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = ASn_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + ASn(i);
    AS2_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = AS2_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + AS2(i);
    u_pert_s2tRMS_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:,:) = u_pert_s2tRMS_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:,:) + u_pert_InNormal_s2tRMS(i,:,:,:);
    uS_pert_s2tRMS_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:,:) = uS_pert_s2tRMS_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:,:) + uS_pert_InNormal_s2tRMS(i,:,:,:);
    u1_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) = u1_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) + u1_InNormal_s2tAvg(i,:,:);
    u1Plus_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) = u1Plus_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) + u1Plus(i,:,:);
    Du1Plus_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) = Du1Plus_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) + Du1Plus(i,:,:);
    uS1_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) = uS1_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) + uS1_InNormal_s2tAvg(i,:,:);
    pPert_s2tRMS_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) = pPert_s2tRMS_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) + pPert_InNormal_s2tRMS(i,:,:);
    nDist_InNormal_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) = nDist_InNormal_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i),:,:) + nDist_InNormal(i,:,:);
    pe_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = pe_stAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + pe_s2tAvg(i);
    ueSmag_s2tAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = ueSmag_s2tAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + ueSmag_s2tAvg(i);
    ueMag_s2tAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = ueMag_s2tAvg_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + ueMag_s2tAvg(i);
    u_tau_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = u_tau_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + u_tau(i);
    yPlus_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = yPlus_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + yPlus(i);
    h_n_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = h_n_tmp(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + h_n(i);
    numPointsZ(mapping1stLayerPoints2Uniq1stLayerPoints(i)) = numPointsZ(mapping1stLayerPoints2Uniq1stLayerPoints(i)) + 1;
end
if min(numPointsZ) < 1; error('Some points were not mapped properly.'); end
delta_s2tAvg = delta_stAvg_tmp ./ numPointsZ;
theta_s2tAvg = theta_stAvg_tmp ./ numPointsZ;
H_s2tAvg = H_stAvg_tmp ./ numPointsZ;
A1 = A1_tmp ./ numPointsZ;
An = An_tmp ./ numPointsZ;
A2 = A2_tmp ./ numPointsZ;
deltaS_s2tAvg = deltaS_stAvg_tmp ./ numPointsZ;
thetaS_s2tAvg = thetaS_stAvg_tmp ./ numPointsZ;
HS_s2tAvg = HS_stAvg_tmp ./ numPointsZ;
AS1 = AS1_tmp ./ numPointsZ;
ASn = ASn_tmp ./ numPointsZ;
AS2 = AS2_tmp ./ numPointsZ;
u_pert_InNormal_s2tRMS = u_pert_s2tRMS_tmp ./ repmat(numPointsZ,[1, maxElemInNormalDir, numPointsNormalDirPerElem, nd]);
uS_pert_InNormal_s2tRMS = uS_pert_s2tRMS_tmp ./ repmat(numPointsZ,[1, maxElemInNormalDir, numPointsNormalDirPerElem, nd]);
u1_InNormal_s2tAvg = u1_stAvg_tmp ./ repmat(numPointsZ, [1, maxElemInNormalDir, numPointsNormalDirPerElem]);
uS1_InNormal_s2tAvg = uS1_stAvg_tmp ./ repmat(numPointsZ, [1, maxElemInNormalDir, numPointsNormalDirPerElem]);
pPert_InNormal_s2tRMS = pPert_s2tRMS_tmp ./ repmat(numPointsZ,[1, maxElemInNormalDir, numPointsNormalDirPerElem]);
nDist_InNormal = nDist_InNormal_tmp ./ repmat(numPointsZ,[1, maxElemInNormalDir, numPointsNormalDirPerElem]);
pe_s2tAvg = pe_stAvg_tmp ./ numPointsZ;
ueSmag_s2tAvg = ueSmag_s2tAvg_tmp ./ numPointsZ;
ueMag_s2tAvg = ueMag_s2tAvg_tmp ./ numPointsZ;
u_tau = u_tau_tmp ./ numPointsZ;
u1Plus = u1Plus_tmp ./ repmat(numPointsZ, [1, maxElemInNormalDir, numPointsNormalDirPerElem]);
Du1Plus = Du1Plus_tmp ./ repmat(numPointsZ, [1, maxElemInNormalDir, numPointsNormalDirPerElem]);
yPlus = yPlus_tmp ./ numPointsZ;
h_n = h_n_tmp ./ numPointsZ;

clear mapping numPointsZ
clear delta_stAvg_tmp theta_stAvg_tmp H_stAvg_tmp A1_tmp An_tmp A2_tmp
clear deltaS_stAvg_tmp thetaS_stAvg_tmp HS_stAvg_tmp AS1_tmp ASn_tmp AS2_tmp
clear nDist_InNormal_tmp pPert_s2tRMS_tmp u1_stAvg_tmp u_pert_s2tRMS_tmp
clear ueSmag_s2tAvg_tmp pe_stAvg_tmp



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get amplification factors of instabilities %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: The leading edge is chosen as the reference location
% TODO: In reality, it should be the point in which instabilities start
% growing, but this might be ok...

N1 = log(A1/min(A1(:)));
N2 = log(A2/min(A2(:)));
NS1 = log(AS1/min(AS1(:)));
NS2 = log(AS2/min(AS2(:)));



%%%%%%%%%%%%%%%%%%%%
%%% Plot results %%%
%%%%%%%%%%%%%%%%%%%%

close all

figure(1)
plot(sFieldLowerSide,delta_s2tAvg(iLowerSide),'.')
hold on
plot(sFieldUpperSide,delta_s2tAvg(iUpperSide),'.')
xlim([-0.1,1.1]);
ylim([0,max(delta_s2tAvg)*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$\delta_1$','Interpreter','latex','FontSize',16)
title('Displacement thickness','Interpreter','latex','FontSize',16)
legend1 = legend('Pressure side','Suction side');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

figure(2)
plot(sFieldLowerSide,theta_s2tAvg(iLowerSide),'.')
hold on
plot(sFieldUpperSide,theta_s2tAvg(iUpperSide),'.')
xlim([-0.1,1.1]);
ylim([0,max(theta_s2tAvg)*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$\theta_1$','Interpreter','latex','FontSize',16)
title('Momentum thickness','Interpreter','latex','FontSize',16)
legend1 = legend('Pressure side','Suction side');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

figure(3)
plot(sFieldLowerSide,H_s2tAvg(iLowerSide),'.')
hold on
plot(sFieldUpperSide,H_s2tAvg(iUpperSide),'.')
xlim([-0.1,1.1]);
ylim([0,max(H_s2tAvg)*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$H_1$','Interpreter','latex','FontSize',16)
title('Shape parameter','Interpreter','latex','FontSize',16)
legend1 = legend('Pressure side','Suction side');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

figure(4)
plot(sFieldLowerSide,N1(iLowerSide),'.')
hold on
plot(sFieldUpperSide,N1(iUpperSide),'.')
legend1 = legend('Pressure side','Suction side');
xlim([-0.1,1.1]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$N_1$','Interpreter','latex','FontSize',16)
title('Amplification factor of streamwise instabilities','Interpreter','latex','FontSize',16)
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

figure(5)
plot(sFieldLowerSide,N2(iLowerSide),'.')
hold on
plot(sFieldUpperSide,N2(iUpperSide),'.')
xlim([-0.1,1.1]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$N_2$','Interpreter','latex','FontSize',16)
title('Amplification factor of cross-flow instabilities','Interpreter','latex','FontSize',16)
legend1 = legend('Pressure side','Suction side');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

u_pert_InNormal_s2tRMS = permute(u_pert_InNormal_s2tRMS, [1,3,2,4]);
u_pert_InNormal_s2tRMS = reshape(u_pert_InNormal_s2tRMS, [], numPointsNormalDirPerElem*maxElemInNormalDir, nd);
uS_pert_InNormal_s2tRMS = permute(uS_pert_InNormal_s2tRMS, [1,3,2,4]);
uS_pert_InNormal_s2tRMS = reshape(uS_pert_InNormal_s2tRMS, [], numPointsNormalDirPerElem*maxElemInNormalDir, nd);
u1_InNormal_s2tAvg = permute(u1_InNormal_s2tAvg, [1,3,2]);
u1_InNormal_s2tAvg = reshape(u1_InNormal_s2tAvg,[], numPointsNormalDirPerElem*maxElemInNormalDir);
uS1_InNormal_s2tAvg = permute(uS1_InNormal_s2tAvg, [1,3,2]);
uS1_InNormal_s2tAvg = reshape(uS1_InNormal_s2tAvg,[], numPointsNormalDirPerElem*maxElemInNormalDir);
u1Plus = permute(u1Plus, [1,3,2]);
u1Plus = reshape(u1Plus,[], numPointsNormalDirPerElem*maxElemInNormalDir);
Du1Plus = permute(Du1Plus, [1,3,2]);
Du1Plus = reshape(Du1Plus,[], numPointsNormalDirPerElem*maxElemInNormalDir);
pPert_InNormal_s2tRMS = permute(pPert_InNormal_s2tRMS, [1,3,2]);
pPert_InNormal_s2tRMS = reshape(pPert_InNormal_s2tRMS, [], numPointsNormalDirPerElem*maxElemInNormalDir);
pe_s2tAvg = pe_s2tAvg(:);

nDist_InNormal = permute(nDist_InNormal, [1,3,2]);
nDist_InNormal = reshape(nDist_InNormal,[], numPointsNormalDirPerElem*maxElemInNormalDir);

inLower = cell(length(iUniqueSfieldLower),1);
inUpper = cell(length(iUniqueSfieldUpper),1);
for i=1:length(iUniqueSfieldLower)
    inLower{i} = find(squeeze(u1_InNormal_s2tAvg(iUniqueSfieldLower(i),:)) ~= 0);
end
for i=1:length(iUniqueSfieldUpper)
    inUpper{i} = find(squeeze(u1_InNormal_s2tAvg(iUniqueSfieldUpper(i),:)) ~= 0);
end

in25PctLower = find(squeeze(u1_InNormal_s2tAvg(iUnique25PctLower,:)) ~= 0);
in50PctLower = find(squeeze(u1_InNormal_s2tAvg(iUnique50PctLower,:)) ~= 0);
in75PctLower = find(squeeze(u1_InNormal_s2tAvg(iUnique75PctLower,:)) ~= 0);
in25PctUpper = find(squeeze(u1_InNormal_s2tAvg(iUnique25PctUpper,:)) ~= 0);
in50PctUpper = find(squeeze(u1_InNormal_s2tAvg(iUnique50PctUpper,:)) ~= 0);
in75PctUpper = find(squeeze(u1_InNormal_s2tAvg(iUnique75PctUpper,:)) ~= 0);

% (s2,t)-RMS of fluctuating streamwise velocity at different BL locations:
% % % % figure(6)
% % % % plot(u_pert_InNormal_s2tRMS(iUnique25PctLower,in25PctLower,1)/ueMag_s2tAvg(iUnique25PctLower), nDist_InNormal(iUnique25PctLower,in25PctLower)/max(nDist_InNormal(iUnique25PctLower,in25PctLower)),'.');
% % % % hold on
% % % % plot(u_pert_InNormal_s2tRMS(iUnique50PctLower,in50PctLower,1)/ueMag_s2tAvg(iUnique50PctLower), nDist_InNormal(iUnique50PctLower,in50PctLower)/max(nDist_InNormal(iUnique50PctLower,in50PctLower)),'.');
% % % % hold on
% % % % plot(u_pert_InNormal_s2tRMS(iUnique75PctLower,in75PctLower,1)/ueMag_s2tAvg(iUnique75PctLower), nDist_InNormal(iUnique75PctLower,in75PctLower)/max(nDist_InNormal(iUnique75PctLower,in75PctLower)),'.');
% % % % hold on
% % % % plot(u_pert_InNormal_s2tRMS(iUnique25PctUpper,in25PctUpper,1)/ueMag_s2tAvg(iUnique25PctUpper), nDist_InNormal(iUnique25PctUpper,in25PctUpper)/max(nDist_InNormal(iUnique25PctUpper,in25PctUpper)),'.');
% % % % hold on
% % % % plot(u_pert_InNormal_s2tRMS(iUnique50PctUpper,in50PctUpper,1)/ueMag_s2tAvg(iUnique50PctUpper), nDist_InNormal(iUnique50PctUpper,in50PctUpper)/max(nDist_InNormal(iUnique50PctUpper,in50PctUpper)),'.');
% % % % hold on
% % % % plot(u_pert_InNormal_s2tRMS(iUnique75PctUpper,in75PctUpper,1)/ueMag_s2tAvg(iUnique75PctUpper), nDist_InNormal(iUnique75PctUpper,in75PctUpper)/max(nDist_InNormal(iUnique75PctUpper,in75PctUpper)),'.');
% % % % xlabel('$\bar{u^{''2}_1} / \bar{u^2_{1,e}}$','Interpreter','latex','FontSize',16);
% % % % ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
% % % % title('Fluctuating streamwise velocity profile at different boundary layer locations','Interpreter','latex','FontSize',16);
% % % % legend1 = legend('Microphone No. 1', 'Microphone No. 2', 'Microphone No. 3', 'Microphone No. 4', 'Microphone No. 5', 'Microphone No. 6');
% % % % set(legend1,'FontSize',14,'Interpreter','latex');
% % % % set(gcf,'color','w');
figure(6)
subplot(1,2,1)
for i=length(iUniqueSfieldLower):-5:1
plot(u_pert_InNormal_s2tRMS(iUniqueSfieldLower(i),inLower{i},1)/ueMag_s2tAvg(iUniqueSfieldLower(i)), nDist_InNormal(iUniqueSfieldLower(i),inLower{i})/max(nDist_InNormal(iUniqueSfieldLower(i),inLower{i})),'.');
hold on
end
xlabel('$\bar{u^{''2}_1} / \bar{u^2_{1,e}}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Fluctuating streamwise velocity profile at different lower BL locations','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

subplot(1,2,2)
for i=1:5:length(iUniqueSfieldUpper)
plot(u_pert_InNormal_s2tRMS(iUniqueSfieldUpper(i),inUpper{i},1)/ueMag_s2tAvg(iUniqueSfieldUpper(i)), nDist_InNormal(iUniqueSfieldUpper(i),inUpper{i})/max(nDist_InNormal(iUniqueSfieldUpper(i),inUpper{i})),'.');
hold on
end
xlabel('$\bar{u^{''2}_1} / \bar{u^2_{1,e}}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Fluctuating streamwise velocity profile at different upper BL locations','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

% Average u1-velocity profile at different BL locations:
figure(7)
plot(u1_InNormal_s2tAvg(iUnique25PctLower,in25PctLower)/u1_InNormal_s2tAvg(iUnique25PctLower,in25PctLower(end)), nDist_InNormal(iUnique25PctLower,in25PctLower)/max(nDist_InNormal(iUnique25PctLower,in25PctLower)),'.');
hold on
plot(u1_InNormal_s2tAvg(iUnique50PctLower,in50PctLower)/u1_InNormal_s2tAvg(iUnique50PctLower,in50PctLower(end)), nDist_InNormal(iUnique50PctLower,in50PctLower)/max(nDist_InNormal(iUnique50PctLower,in50PctLower)),'.');
hold on
plot(u1_InNormal_s2tAvg(iUnique75PctLower,in75PctLower)/u1_InNormal_s2tAvg(iUnique75PctLower,in75PctLower(end)), nDist_InNormal(iUnique75PctLower,in75PctLower)/max(nDist_InNormal(iUnique75PctLower,in75PctLower)),'.');
hold on
plot(u1_InNormal_s2tAvg(iUnique25PctUpper,in25PctUpper)/u1_InNormal_s2tAvg(iUnique25PctUpper,in25PctUpper(end)), nDist_InNormal(iUnique25PctUpper,in25PctUpper)/max(nDist_InNormal(iUnique25PctUpper,in25PctUpper)),'.');
hold on
plot(u1_InNormal_s2tAvg(iUnique50PctUpper,in50PctUpper)/u1_InNormal_s2tAvg(iUnique50PctUpper,in50PctUpper(end)), nDist_InNormal(iUnique50PctUpper,in50PctUpper)/max(nDist_InNormal(iUnique50PctUpper,in50PctUpper)),'.');
hold on
plot(u1_InNormal_s2tAvg(iUnique75PctUpper,in75PctUpper)/u1_InNormal_s2tAvg(iUnique75PctUpper,in75PctUpper(end)), nDist_InNormal(iUnique75PctUpper,in75PctUpper)/max(nDist_InNormal(iUnique75PctUpper,in75PctUpper)),'.');
xlim([-0.1,1.1]); ylim([0,1]);
xlabel('$\bar{u}_1 / \bar{u}_{1,e}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Average streamwise-velocity profile at different boundary layer locations','Interpreter','latex','FontSize',16);
legend1 = legend('Microphone No. 1', 'Microphone No. 2', 'Microphone No. 3', 'Microphone No. 4', 'Microphone No. 5', 'Microphone No. 6');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

% (s2,t)-RMS of fluctuating pressure at different BL locations:
% % % % % % figure(8)
% % % % % % plot(pPert_InNormal_s2tRMS(iUnique25PctLower,in25PctLower) / pe_s2tAvg(iUnique25PctLower), nDist_InNormal(iUnique25PctLower,in25PctLower)/nDist_InNormal(iUnique25PctLower,in25PctLower(end)),'.');
% % % % % % hold on
% % % % % % plot(pPert_InNormal_s2tRMS(iUnique50PctLower,in50PctLower) / pe_s2tAvg(iUnique50PctLower), nDist_InNormal(iUnique50PctLower,in50PctLower)/nDist_InNormal(iUnique50PctLower,in50PctLower(end)),'.');
% % % % % % hold on
% % % % % % plot(pPert_InNormal_s2tRMS(iUnique75PctLower,in75PctLower) / pe_s2tAvg(iUnique75PctLower), nDist_InNormal(iUnique75PctLower,in75PctLower)/nDist_InNormal(iUnique75PctLower,in75PctLower(end)),'.');
% % % % % % hold on
% % % % % % plot(pPert_InNormal_s2tRMS(iUnique25PctUpper,in25PctUpper) / pe_s2tAvg(iUnique25PctUpper), nDist_InNormal(iUnique25PctUpper,in25PctUpper)/nDist_InNormal(iUnique25PctUpper,in25PctUpper(end)),'.');
% % % % % % hold on
% % % % % % plot(pPert_InNormal_s2tRMS(iUnique50PctUpper,in50PctUpper) / pe_s2tAvg(iUnique50PctUpper), nDist_InNormal(iUnique50PctUpper,in50PctUpper)/nDist_InNormal(iUnique50PctUpper,in50PctUpper(end)),'.');
% % % % % % hold on
% % % % % % plot(pPert_InNormal_s2tRMS(iUnique75PctUpper,in75PctUpper) / pe_s2tAvg(iUnique75PctUpper), nDist_InNormal(iUnique75PctUpper,in75PctUpper)/nDist_InNormal(iUnique75PctUpper,in75PctUpper(end)),'.');
% % % % % % xlabel('$p^{''}_{RMS} / p_e$','Interpreter','latex','FontSize',16);
% % % % % % ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
% % % % % % title('Fluctuating pressure profile at different boundary layer locations','Interpreter','latex','FontSize',16);
% % % % % % legend1 = legend('Microphone No. 1', 'Microphone No. 2', 'Microphone No. 3', 'Microphone No. 4', 'Microphone No. 5', 'Microphone No. 6');
% % % % % % set(legend1,'FontSize',14,'Interpreter','latex');
% % % % % % set(gcf,'color','w');
figure(8)
subplot(1,2,1)
for i=length(iUniqueSfieldLower):-5:1
plot(pPert_InNormal_s2tRMS(iUniqueSfieldLower(i),inLower{i}) / pe_s2tAvg(iUniqueSfieldLower(i)), nDist_InNormal(iUniqueSfieldLower(i),inLower{i})/max(nDist_InNormal(iUniqueSfieldLower(i),inLower{i})),'.');
hold on
end
xlabel('$\bar{u^{''2}_1} / \bar{u^2_{1,e}}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Fluctuating pressure profile at different lower BL locations','Interpreter','latex','FontSize',16);
set(gcf,'color','w');

subplot(1,2,2)
for i=1:5:length(iUniqueSfieldUpper)
plot(pPert_InNormal_s2tRMS(iUniqueSfieldUpper(i),inUpper{i}) / pe_s2tAvg(iUniqueSfieldUpper(i)), nDist_InNormal(iUniqueSfieldUpper(i),inUpper{i})/max(nDist_InNormal(iUniqueSfieldUpper(i),inUpper{i})),'.');
hold on
end
xlabel('$\bar{u^{''2}_1} / \bar{u^2_{1,e}}$','Interpreter','latex','FontSize',16);
ylabel('$n / n_e$','Interpreter','latex','FontSize',16);
title('Fluctuating pressure profile at different upper BL locations','Interpreter','latex','FontSize',16);
set(gcf,'color','w');



% % Fourier transform of pressure and turbulent kinetic energy at different BL locations:
% figure(9)
% subplot(1,2,1)
% for i=1:numMicrophones
%     hold on
%     plot(frequencies,FFTpAtMicrophones(:,i),'.');
% end
% xlabel('Frequency $f$ (Hz)','Interpreter','latex','FontSize',16);
% ylabel('Amplitude of FFT($p$)','Interpreter','latex','FontSize',16);
% title('FFT of pressure at microphone locations','Interpreter','latex','FontSize',16);
% legend1 = legend('Microphone No. 1', 'Microphone No. 2', 'Microphone No. 3', 'Microphone No. 4', 'Microphone No. 5', 'Microphone No. 6');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');
% 
% subplot(1,2,2)
% for i=1:numMicrophones
%     hold on
%     plot(frequencies,FFTkAtMicrophones(:,i),'.');
% end
% xlabel('Frequency $f$ (Hz)','Interpreter','latex','FontSize',16);
% ylabel('Amplitude of FFT($k$)','Interpreter','latex','FontSize',16);
% title('FFT of turbulent kinetic energy at microphone locations','Interpreter','latex','FontSize',16);
% legend1 = legend('Microphone No. 1', 'Microphone No. 2', 'Microphone No. 3', 'Microphone No. 4', 'Microphone No. 5', 'Microphone No. 6');
% set(legend1,'FontSize',14,'Interpreter','latex');
% set(gcf,'color','w');



% Plot y+ and h_n:
h_n_lowerSide = h_n(iLowerSide);
h_n_upperSide = h_n(iUpperSide);
yPlus_lowerSide = yPlus(iLowerSide);
yPlus_upperSide = yPlus(iUpperSide);

figure(10)
subplot(1,2,1)
plot(sFieldLowerSide,h_n_lowerSide,'.')
hold on
plot(sFieldUpperSide,h_n_upperSide,'.')
xlim([-0.1,1.1]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$h_n$','Interpreter','latex','FontSize',16)
title('Normal element size near aifoil','Interpreter','latex','FontSize',16)
legend1 = legend('Pressure side','Suction side');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

subplot(1,2,2)
plot(sFieldLowerSide,yPlus_lowerSide,'.')
hold on
plot(sFieldUpperSide,yPlus_upperSide,'.')
xlim([-0.1,1.1]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$y^{+}$','Interpreter','latex','FontSize',16)
title('$y^{+}$ for first cell','Interpreter','latex','FontSize',16)
legend1 = legend('Pressure side','Suction side');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');



% [gaussPoints1d_sWise,gaussWeights1d_sWise] = gaussquad(pgauss_sWise,1,0);
clear m l
H_s2tAvg_min = 1.5;
H_s2tAvg_max = 20;
H_s2tAvg_reg = min(max(H_s2tAvg,H_s2tAvg_min),H_s2tAvg_max);

% Ref.: M. Drela and M. Giles. Viscous-inviscid analysis of transonic and low Reynolds number airfoils. AIAA Journal, 25(10):1347–1355, 1987
dNdRe_th_v1 = 0.01 * sqrt((2.4*H_s2tAvg_reg - 3.7 + 2.5 * tanh(1.5*H_s2tAvg_reg - 4.65)).^2 + 0.25);
l = (6.54*H_s2tAvg_reg - 14.07) ./ H_s2tAvg_reg.^2;
m = (0.058*((H_s2tAvg_reg-4).^2./(H_s2tAvg_reg-1)) - 0.068) .* (1./l);
dRe_th_dx_v1 = ((m+1)/2) .* l ./ theta_s2tAvg;
Re_theta_cr_lowerSide_v1 = 10.^(((1.415./(H_s2tAvg_reg(iLowerSide)-1)) - 0.489).*tanh(20./(H_s2tAvg_reg(iLowerSide)-1)-12.9) + 3.295./(H_s2tAvg_reg(iLowerSide)-1)+0.44);
Re_theta_cr_upperSide_v1 = 10.^(((1.415./(H_s2tAvg_reg(iUpperSide)-1)) - 0.489).*tanh(20./(H_s2tAvg_reg(iUpperSide)-1)-12.9) + 3.295./(H_s2tAvg_reg(iUpperSide)-1)+0.44);

% Ref.: M. Drela notes for 16.13 course
dNdRe_th_v2 = 0.028*(H_s2tAvg_reg-1)-0.0345*exp(-(3.87./(H_s2tAvg_reg-1)-2.52).^2);
dRe_th_dx_v2 = (-0.05+2.7./(H_s2tAvg_reg-1)-5.5./(H_s2tAvg_reg-1).^2+3.0./(H_s2tAvg_reg-1).^3+0.1*exp(-20./(H_s2tAvg_reg-1))) ./ theta_s2tAvg;
Re_theta_cr_lowerSide_v2 = 10.^(2.492./(H_s2tAvg_reg(iLowerSide)-1).^0.43 + 0.7*(tanh(14./(H_s2tAvg_reg(iLowerSide)-1)-9.24)+1));
Re_theta_cr_upperSide_v2 = 10.^(2.492./(H_s2tAvg_reg(iUpperSide)-1).^0.43 + 0.7*(tanh(14./(H_s2tAvg_reg(iUpperSide)-1)-9.24)+1));

if eNformulation == 1       % M. Drela and M. Giles. Viscous-inviscid analysis of transonic and low Reynolds number airfoils. AIAA Journal, 25(10):1347–1355, 1987
    dNdRe_th = dNdRe_th_v1;
    dRe_th_dx = dRe_th_dx_v1;
    Re_theta_cr_lowerSide = Re_theta_cr_lowerSide_v1;
    Re_theta_cr_upperSide = Re_theta_cr_upperSide_v1;
elseif eNformulation == 2   % M. Drela notes for 16.13 course
    dNdRe_th = dNdRe_th_v2;
    dRe_th_dx = dRe_th_dx_v2;
    Re_theta_cr_lowerSide = Re_theta_cr_lowerSide_v2;
    Re_theta_cr_upperSide = Re_theta_cr_upperSide_v2;
else
    error('eN formulation flag has invalid value.');
end

integrand = dNdRe_th .* dRe_th_dx;

Re_theta_lowerSide = Re_inf * thetaS_s2tAvg(iLowerSide) .* ueSmag_s2tAvg(iLowerSide);
Re_theta_upperSide = Re_inf * thetaS_s2tAvg(iUpperSide) .* ueSmag_s2tAvg(iUpperSide);


% Get unstable portion of BLs
iLowerSideUnstable_aux = find(Re_theta_lowerSide(:) > Re_theta_cr_lowerSide(:),10,'last');
iLowerSideUnstable_aux = flip(iLowerSideUnstable_aux);
for i=1:10
    if all(Re_theta_lowerSide((iLowerSideUnstable_aux(i)-10):iLowerSideUnstable_aux(i)) > Re_theta_cr_lowerSide((iLowerSideUnstable_aux(i)-10):iLowerSideUnstable_aux(i)))
        iLowerSideUnstable_aux = iLowerSideUnstable_aux(i);
        break
    end
end
if i==10; error('Starting point of unstable BL was not properly detected.'); end
iLowerSideStable = (iLowerSideUnstable_aux+1):length(iLowerSide);
iLowerSideUnstable = 1:iLowerSideUnstable_aux;
iUpperSideUnstable_aux = find(Re_theta_upperSide(:) > Re_theta_cr_upperSide(:),10,'first');
for i=1:10
    if all(Re_theta_upperSide(iUpperSideUnstable_aux(i):(iUpperSideUnstable_aux(i)+10)) > Re_theta_cr_upperSide(iUpperSideUnstable_aux(i):(iUpperSideUnstable_aux(i)+10)))
        iUpperSideUnstable_aux = iUpperSideUnstable_aux(i);
        break
    end
end
if i==10; error('Starting point of unstable BL was not properly detected.'); end
iUpperSideStable = length(iLowerSide)-1 + (1:(iUpperSideUnstable_aux-1));
iUpperSideUnstable = length(iLowerSide)-1 + (iUpperSideUnstable_aux:length(iUpperSide));
clear iLowerSideUnstable_aux iUpperSideUnstable_aux


% Get turbulent portion of BLs
N1_lowerSide = N1(iLowerSide);
N1_upperSide = N1(iUpperSide);
iLowerSideTurbulent_aux = find(N1_lowerSide > N1_cr_lowerSide, 10, 'last');
iLowerSideTurbulent_aux = flip(iLowerSideTurbulent_aux);
for i=1:10
    if all(N1_lowerSide((iLowerSideTurbulent_aux(i)-5):iLowerSideTurbulent_aux(i)) > N1_cr_lowerSide)
        iLowerSideTurbulent_aux = iLowerSideTurbulent_aux(i);
        break
    end
end
if i==10; error('Starting point of turbulent lower BL was not properly detected.'); end
iLowerSideTurbulent = 1:iLowerSideTurbulent_aux;
iUpperSideTurbulent_aux = find(N1_upperSide > N1_cr_upperSide,10,'first');
for i=1:10
    if all(N1_upperSide(iUpperSideTurbulent_aux(i):(iUpperSideTurbulent_aux(i)+5)) > N1_cr_upperSide)
        iUpperSideTurbulent_aux = iUpperSideTurbulent_aux(i);
        break
    end
end
if i==10; error('Starting point of turbulent upper BL was not properly detected.'); end
iUpperSideTurbulent = length(iLowerSide)-1 + (iUpperSideTurbulent_aux:length(iUpperSide));
clear iLowerSideTurbulent_aux iUpperSideTurbulent_aux

% sFieldLowerSideUnstable = sField(iLowerSideUnstable);
% sFieldUpperSideUnstable = sField(iUpperSideUnstable);
% dNdRe_th_lowerSideUnstable = dNdRe_th(iLowerSideUnstable);
% dNdRe_th_upperSideUnstable = dNdRe_th(iUpperSideUnstable);
integrand_lowerSideUnstable = integrand(iLowerSideUnstable);
integrand_upperSideUnstable = integrand(iUpperSideUnstable);

if numElemOnAirfoilPerExtrLayer*numPointsPerElem_x ~= size(xyCoordUnique,1); error('Error.'); end

% weights1d_s1Wise = (1/numPointsPerElem_x) * ones(numPointsPerElem_x,1);
% weights1d_s1Wise = repmat(weights1d_s1Wise,[numElemOnAirfoilPerExtrLayer,1]);
% weights1d_s1Wise_lowerSideUnstable = weights1d_s1Wise(iLowerSideUnstable);
% weights1d_s1Wise_upperSideUnstable = weights1d_s1Wise(iUpperSideUnstable);

ds_lowerSideUnstable_tmp = sqrt(sum((xyCoordUnique(iLowerSideUnstable(2:end),1:2) - xyCoordUnique(iLowerSideUnstable(1:end-1),1:2)).^2,2));       % TODO: Compute exact length
ds_lowerSideUnstable_tmp = ds_lowerSideUnstable_tmp(:);
ds_lowerSideUnstable = 0.5*ds_lowerSideUnstable_tmp(1);
ds_lowerSideUnstable = [ds_lowerSideUnstable; 0.5*ds_lowerSideUnstable_tmp(1:end-1) + 0.5*ds_lowerSideUnstable_tmp(2:end); 0.5*ds_lowerSideUnstable_tmp(end)];

ds_upperSideUnstable_tmp = sqrt(sum((xyCoordUnique(iUpperSideUnstable(2:end),1:2) - xyCoordUnique(iUpperSideUnstable(1:end-1),1:2)).^2,2));       % TODO: Compute exact length
ds_upperSideUnstable_tmp = ds_upperSideUnstable_tmp(:);
ds_upperSideUnstable = 0.5*ds_upperSideUnstable_tmp(1);
ds_upperSideUnstable = [ds_upperSideUnstable; 0.5*ds_upperSideUnstable_tmp(1:end-1) + 0.5*ds_upperSideUnstable_tmp(2:end); 0.5*ds_upperSideUnstable_tmp(end)];

% Integration from leading edge:
% TODO: Implement a more accurate numerical quadrature rule.
N_envelopeMethod_lowerSideUnstable = flip(cumsum(flip(integrand_lowerSideUnstable.*ds_lowerSideUnstable)));
N_envelopeMethod_upperSideUnstable = cumsum(integrand_upperSideUnstable.*ds_upperSideUnstable);

N_envelopeMethod_lowerSide = [N_envelopeMethod_lowerSideUnstable(:); zeros(length(iLowerSideStable),1)];
N_envelopeMethod_upperSide = [zeros(length(iUpperSideStable),1); N_envelopeMethod_upperSideUnstable(:)];

figure(11)
plot(sFieldLowerSide, N1(iLowerSide),'.')
hold on
plot(sFieldUpperSide, N1(iUpperSide),'.')
hold on
plot(sFieldLowerSide, N_envelopeMethod_lowerSide,'.');
hold on
plot(sFieldUpperSide, N_envelopeMethod_upperSide,'.');
xlim([-0.1,1.1]);
% ylim([0,max(max(N_envelopeMethod_lowerSide(:)),max(N_envelopeMethod_upperSide(:)))*1.15]);
xlabel('$x/c$','Interpreter','latex','FontSize',16)
ylabel('$N_1$','Interpreter','latex','FontSize',16)
title('Amplification factor of streamwise instabilities','Interpreter','latex','FontSize',16)
legend1 = legend('LES - Pressure side','LES - Suction side','Envelope method - Pressure side','Envelope method - Suction side');
set(legend1,'FontSize',14,'Interpreter','latex');
set(gcf,'color','w');

clear N_envelopeMethod_lowerSideUnstable N_envelopeMethod_upperSideUnstable
clear ds_upperSideUnstable ds_upperSideUnstable_tmp
clear ds_lowerSideUnstable ds_lowerSideUnstable_tmp



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute resolved flow features %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elemNextToAirfoil1stExtrLayer = zeros(numFacesOnAirfoilFirstExtrLayer,1);
for i=1:numFacesOnAirfoilFirstExtrLayer
    elemNextToAirfoil1stExtrLayer(i) = cell2mat(elemInNormalDir{facesOnAirfoilFirstExtrLayer_aux(i),1}(1));
end
hFieldElemNextToAirfoil1stExtrLayer = hField(:,elemNextToAirfoil1stExtrLayer);
hFieldElemNextToAirfoil1stExtrLayer_spanwise = squeeze(hFieldElemNextToAirfoil1stExtrLayer(1,:));
hFieldElemNextToAirfoil1stExtrLayer_normal = squeeze(hFieldElemNextToAirfoil1stExtrLayer(2,:));
hFieldElemNextToAirfoil1stExtrLayer_streamwise = squeeze(hFieldElemNextToAirfoil1stExtrLayer(3,:));
maxResolvablePhaseSpeedTSwaves_totalBL = (min(hFieldElemNextToAirfoil1stExtrLayer_streamwise)/porder) / (dt / torder);
maxResolvablePhaseSpeedTSwaves_unstableBL = (min(hFieldElemNextToAirfoil1stExtrLayer_streamwise([iLowerSideUnstable(:); iUpperSideUnstable(:)]))/porder) / (dt / torder);
disp(' ');
disp(['Minimum resolvable wavelength for TS waves is (Nyquist-Shannon theorem): ', num2str(2*max(hFieldElemNextToAirfoil1stExtrLayer_streamwise([iLowerSideUnstable(:); iUpperSideUnstable(:)]))/porder)]);
disp(['Maximum resolvable phase speed of TS waves is: ', num2str(maxResolvablePhaseSpeedTSwaves_unstableBL)]);

minSpaceResolutionIsotropicTurbulence = max(reshape(hFieldElemNextToAirfoil1stExtrLayer(:,[iLowerSideTurbulent(:); iUpperSideTurbulent(:)]),[],1))/porder;
minTimeResolutionTurbulence = dt/torder;
disp(' ');
disp(['Minimum resolvable wavelength of isotropic turbulence is (Nyquist-Shannon theorem): ', num2str(2*minSpaceResolutionIsotropicTurbulence)]);
disp(['Minimum resolvable time period of turbulent eddies is (Nyquist-Shannon theorem): ', num2str(2*minTimeResolutionTurbulence)]);
if elemtype ~= 1; warning('All indicators of resolver flow features are only accurate for hex elements'); end


%%%%%%%%%%%%%%%%%%
%%% Clear data %%%
%%%%%%%%%%%%%%%%%%

clear r ru rv rw r_x rv_x rw_x r_y ru_y rw_y r_z ru_z rv_z
clear u v w v_x w_x u_y w_y u_z v_z
clear vort_x Dvort_xDn vortCrossN vortMag DvortDnMag
clear porder pgauss elemtype plocvl plocfc npf nfe nd perm dgnodes nc npv ne nf nch pMesh tMesh fMesh t2fMesh
clear shapg UtoLineIntegralU elemInNormalDir lenInElem numElemInNormalDir



%%%%%%%%%%%%%%%%%%%%
%%% Save results %%%
%%%%%%%%%%%%%%%%%%%%

end


function xi = mapBack(dgnodes, x, porder, elemtype, plocvl, npv, nd, Alocvl_inv, pp_vl, dpp_vl)

maxDelta = 0.1;
margin = 1e-1;
tol = 1e-8; %1e-9;
maxNumLineSearch = 6; %20;
maxIter = 15; %50;

dgnodes = squeeze(dgnodes);

[numPoints,~] = size(x);

if nd ~= size(x,2); error('x has wrong size at input.'); end

% Get initial guess:
iInitialGuess = zeros(numPoints,1);
for i=1:numPoints
    initialGuess_aux = repmat(x(i,:)',[1,npv]) - dgnodes';
    [~,iInitialGuess(i)] = min(sum(initialGuess_aux.^2,1));
end
initialGuess = plocvl(iInitialGuess,:);      % [numPoints,nd]

xi = zeros(numPoints,nd);
for i=1:numPoints
    % Newton iteration:
    iter = 0;
    while true
        iter = iter + 1;
        alpha = 1.0;

        if iter == 1
            xi_n = initialGuess(i,:);
            xi_n = xi_n(:);
            shap_n = mkshape(porder, plocvl, reshape(xi_n,[1,nd]), elemtype, Alocvl_inv, pp_vl, dpp_vl);
            shap_n = squeeze(shap_n(:,1,:));
            f_n = x(i,:)' - dgnodes'*shap_n(:,1);
        else
            shap_n = shap_nPlus1;
            f_n = f_nPlus1;
            xi_n = xi_nPlus1;
        end
        J_n = - dgnodes'*shap_n(:,2:end);

        delta = J_n\f_n;
        delta = delta(:);
        if norm(delta) > maxDelta; delta = maxDelta * delta / norm(delta); end
        xi_nPlus1 = xi_n - alpha*delta;
        shap_nPlus1 = mkshape(porder, plocvl, reshape(xi_nPlus1,[1,nd]), elemtype, Alocvl_inv, pp_vl, dpp_vl);
        shap_nPlus1 = squeeze(shap_nPlus1(:,1,:));
        f_nPlus1 = x(i,:)' - dgnodes'*shap_nPlus1(:,1);
        isInElement = getIsInElement(xi_nPlus1, elemtype, margin);
        
        % Line search
        numLineSearch = 0;
        while (norm(f_nPlus1) > norm(f_n) &&  norm(f_nPlus1) >= tol) || isInElement == 0 % IT COULD HAPPEN THAT THE INITIAL GUESS IS ON THE BOUNDARY AND THE MAXIMUM DESCEND DIRECTION POINTS OUTWARDS!!
            numLineSearch = numLineSearch + 1;
            alpha = alpha / 2;
            xi_nPlus1 = xi_n - alpha*delta;
            shap_nPlus1 = mkshape(porder, plocvl, reshape(xi_nPlus1,[1,nd]), elemtype, Alocvl_inv, pp_vl, dpp_vl);
            shap_nPlus1 = squeeze(shap_nPlus1(:,1,:));
            f_nPlus1 = x(i,:)' - dgnodes'*shap_nPlus1(:,1);
            isInElement = getIsInElement(xi_nPlus1, elemtype, margin);
            if numLineSearch > maxNumLineSearch;
                warning('Too many line searchs.');
                break;
            end
        end

        if norm(f_nPlus1) < tol
            break;
        end
        
        if iter >= maxIter;
            warning(['Too many Newton iterations. Final residual: ', num2str(norm(f_nPlus1))]);
            break;
        end
    end

    xi(i,:) = xi_nPlus1;
end

end


function isInElement = getIsInElement(xi, elemtype, margin)

if elemtype == 0
    if min(xi) < -margin || sum(xi) > 1+margin
        isInElement = 0;
    else
        isInElement = 1;
    end
elseif elemtype == 1
    if min(xi) < -margin || max(xi) > 1+margin
        isInElement = 0;
    else
        isInElement = 1;
    end
end

end


function [nextFace_local,dist,finalPoint] = getNextFace(normal,initialPoint,currentLocalFace,p,porder,plocfc,dgnodes,perm,elemtype,npf, Alocfc_inv, pp_fc, dpp_fc)

tolPlane = 1e-9;

% Make column vectors:
normal = normal(:);
initialPoint = initialPoint(:);
normalExtended = repmat(normal,[1,npf]);
initialPointExtended = repmat(initialPoint,[1,npf]);

% Get sizes:
nv = size(p,1);     % Number of vertices in element
nd = size(p,2);

if nd == 3
    if nv == 4          % Tetra
        numFaces = 4;
        ind = [[2,3,4];[1,4,3];[1,2,4];[1,3,2]];
    elseif nv == 8      % Hexa
        numFaces = 6;
        ind = [[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]];
    else
        error('Element type not implemented.');
    end
else
    error('Number of dimensions not implemented.');
end

p1 = zeros(nd,numFaces);
p2 = zeros(nd,numFaces);
p3 = zeros(nd,numFaces);
t1 = zeros(nd,numFaces);
t2 = zeros(nd,numFaces);
u = zeros(numFaces,1);
v = zeros(numFaces,1);
t = zeros(numFaces,1);
intersectionPoint = zeros(nd,numFaces);
for i=1:numFaces
    p1(:,i) = p(ind(i,1),:);
    p2(:,i) = p(ind(i,2),:);
    p3(:,i) = p(ind(i,3),:);
    t1(:,i) = p2(:,i) - p1(:,i);
    t2(:,i) = p3(:,i) - p1(:,i);
    
    n_i = cross(t1(:,i),t2(:,i));
    n_i = repmat(n_i(:)'/norm(n_i),[npf,1]);
    
    dgnodesFace = dgnodes(perm(:,i),1:3);
    
    aux = abs(dot( n_i , dgnodesFace-repmat(p1(:,i)',[npf,1]) ,2));
    maxDist = max(aux);
    
    if maxDist < tolPlane       % Linear face
        sol = [t1(:,i), t2(:,i), normal]\(p1(:,i)-initialPoint);
        
        u(i) = sol(1);
        v(i) = sol(2);
        t(i) = sol(3);
        
        finalPoint_v1 = initialPoint + t(i)*normal;
        finalPoint_v2 = p1(:,i) - u(i) * t1(:,i) - v(i) * t2(:,i);
        
        if norm(finalPoint_v1-finalPoint_v2)/max(abs(t(i)),1) > 1.0e-8; error('Computation of line-plane intersection went wrong.'); end
        
        intersectionPoint(:,i) = finalPoint_v1(:);
    else                        % Curved face
        % Get initial guess:
        aux = dot(dgnodesFace' - initialPointExtended, normalExtended, 1);
        aux2 = repmat(initialPoint,[1,npf]) + normal * aux;
        [~, iMinDist] = min(sqrt(sum((dgnodesFace' - aux2).^2,1)));
        initialGuess(1) = aux(iMinDist);
        initialGuess(2:3) = plocfc(iMinDist,:);
        initialGuess = initialGuess(:);
        
        % Compute intersection:
        [t(i),intersectionPoint(:,i),noIntersections] = getIntersectionWithCurvedFace(initialPoint, normal, dgnodesFace, initialGuess, porder, plocfc, elemtype, Alocfc_inv, pp_fc, dpp_fc);
    end
end

if abs(t(currentLocalFace)) > 1.0e-8;
    error('Initial point does not lie on initial face.');
end

% Remove distance to current face
t = [t(1:currentLocalFace-1);t(currentLocalFace+1:end)];

% "Remove" faces with negative distance
toChange = t < 0;
t(toChange) = Inf;

[dist,iFace] = min(t);

if dist == Inf; error('No intersection in the positive normal direction was detected.'); end

if iFace >= currentLocalFace; iFace = iFace + 1; end
nextFace_local = iFace;
finalPoint = intersectionPoint(:,nextFace_local);

end


function [dist,intersectionPoint,noIntersections] = getIntersectionWithCurvedFace(initialPoint, normal, dgnodesFace, initialGuess, porder, plocfc, elemtype, Alocfc_inv, pp_fc, dpp_fc)

maxDelta = 0.1;
margin = 1.0e-1;
tol = 1e-8; %1e-9;
maxNumLineSearch = 6; %20;
maxIter = 15; %50;

normal = normal(:);

% Newton iteration:
iter = 0;
while true
    iter = iter + 1;
    alpha = 1.0;
    
    if iter == 1
        xi_n = initialGuess(:);
        shap_xn = mkshape(porder,plocfc,reshape(xi_n(2:3),[1,2]),elemtype,Alocfc_inv,pp_fc,dpp_fc);
        shap_xn = squeeze(shap_xn(:,1,:));
        f_n = initialPoint + xi_n(1)*normal - dgnodesFace'*shap_xn(:,1);
    else
        shap_xn = shap_nPlus1;
        f_n = f_nPlus1;
        xi_n = xi_nPlus1;
    end
    J_n = [normal, -dgnodesFace'*shap_xn(:,2:3)]; 
    
    delta = J_n\f_n;
    delta = delta(:);
    if norm(delta(2:3)) > maxDelta; delta = maxDelta * delta / norm(delta(2:3)); end
    xi_nPlus1 = xi_n - alpha*delta;
    
    shap_nPlus1 = mkshape(porder, plocfc, reshape(xi_nPlus1(2:3),[1,2]), elemtype, Alocfc_inv, pp_fc, dpp_fc);
    shap_nPlus1 = squeeze(shap_nPlus1(:,1,:));
    f_nPlus1 = initialPoint + xi_nPlus1(1)*normal - dgnodesFace'*shap_nPlus1(:,1);
    isInElement = getIsInElement(xi_nPlus1, elemtype, margin);
    
    % Line search
    numLineSearch = 0;
    while (norm(f_nPlus1) > norm(f_n) &&  norm(f_nPlus1) >= tol) || isInElement == 0 % IT COULD HAPPEN THAT THE INITIAL GUESS IS ON THE BOUNDARY AND THE MAXIMUM DESCEND DIRECTION POINTS OUTWARDS!!
        numLineSearch = numLineSearch + 1;
        alpha = alpha / 2;
        xi_nPlus1 = xi_n - alpha*delta;
        shap_nPlus1 = mkshape(porder, plocfc, reshape(xi_nPlus1(2:3),[1,2]), elemtype, Alocfc_inv, pp_fc, dpp_fc);
        shap_nPlus1 = squeeze(shap_nPlus1(:,1,:));
        f_nPlus1 = initialPoint + xi_nPlus1(1)*normal - dgnodesFace'*shap_nPlus1(:,1);
        isInElement = getIsInElement(xi_nPlus1, elemtype, margin);
        if numLineSearch > maxNumLineSearch;
            break;
        end
    end

    if norm(f_nPlus1) < tol
        noIntersections = 1;
        break;
    end

    if iter >= maxIter % || norm(delta) < 1.0e-8
        noIntersections = 0;
        intersectionPoint = [0; 0; 0];
        dist = Inf;
        break;
    end
end

if noIntersections == 1
    dist = xi_nPlus1(1);
    
    intersectionPoint_v1 = initialPoint + dist * normal;
    intersectionPoint_v2 = dgnodesFace'*shap_nPlus1(:,1);
    if norm(intersectionPoint_v1-intersectionPoint_v2) > tol; error('Computation of intersection between line and polynomial surface went wrong.'); end

    intersectionPoint = intersectionPoint_v1;
end

end


function [coordPointsMaster2d, coordPointsMaster3d] = getCoordPointsMaster(mesh, indexFace, numPoints)

nd = mesh.nd;
elemtype = mesh.elemtype;

if nd ~= 3; error('getCoordPointsMaster only implemented for 3D meshes.'); end

% if elemtype == 0        % Tetra
%     ind = [[2,3,4];[1,4,3];[1,2,4];[1,3,2]];
%     p = [[0,0,0];[1,0,0];[0,1,0];[0,0,1]];
% elseif elemtype == 1    % Hexa
%     ind = [[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]];
%     p = [[0,0,0];[1,0,0];[0,1,0];[1,1,0];[0,0,1];[1,0,1];[0,1,1];[1,1,1]];
% else
%     error('elemtype not implemented.');
% end

% pFace = p(ind(indexFace,:),:);

% positions = (0:numPoints-1) / (numPoints-1);
% coordPointsMaster2d = repmat([0;0],[1,numPoints]) + [1;0] * positions;
% coordPointsMaster2d = repmat(coordPointsMaster2d,[1,numPoints]) + [0; 1] * kron(positions,ones(1,numPoints));

if elemtype == 0
    positions = (0.5:(2*numPoints-0.5))/(2*numPoints);
    coordPointsMaster2d(1,:) = kron(ones(1,2*numPoints),positions);
    coordPointsMaster2d(2,:) = kron(positions,ones(1,2*numPoints));
    pointsToKeep = rot90(tril(ones(2*numPoints,2*numPoints),-1),-1);
    pointsToKeep = pointsToKeep(:);
    pointsToKeep = find(pointsToKeep == 1);
    coordPointsMaster2d = coordPointsMaster2d(:,pointsToKeep);
elseif elemtype == 1
    positions = (0.5:(numPoints-0.5))/numPoints;
    coordPointsMaster2d(1,:) = kron(ones(1,numPoints),positions);
    coordPointsMaster2d(2,:) = kron(positions,ones(1,numPoints));

    % coordPointsMaster3d = repmat(pFace(1,:)',[1,numPoints]) + positions * (pFace(2,:)-pFace(1,:)');
    % coordPointsMaster3d = repmat(coordPointsMaster3d,[1,numPoints]) + kron(positions,ones(numPoints,1)) * (pFace(3,:)-pFace(1,:)');

    coordPointsMaster2d = reshape(coordPointsMaster2d,nd-1,[])';      % [numPoints,nd-1]
end
coordPointsMaster2d = reshape(coordPointsMaster2d,[],2);

coordPointsMaster3d = [];

end


function s = get_s_coord(LEcoord,alpha,x)

[numPoints, nd] = size(x);

if nd ~= 3; error('Wrong number of dimensions detected.'); end

s = zeros(numPoints, 1);

LEcoord = [LEcoord(:); 0];
sVector = [cos(alpha); sin(alpha); 0];
for i=1:numPoints
    s(i) = dot(x(i,:)' - LEcoord, sVector);

    if s(i) < 0; warning('Point on airfoil surface with s < 0 was detected.'); end
end

end


function masterCentroid = getMasterCentroid(elemtype, nd)

if elemtype == 0
    masterCentroid = (1/3) * ones(1,nd);
elseif elemtype == 1
    masterCentroid = 0.5 * ones(1,nd);
end

end


function nfs = mkshape(porder,plocal,pts,elemtype,Ainv,pp,dpp)

if nargin < 5
    if elemtype==0 % simplex elements
        A = koornwinder(plocal,porder,pp,dpp);               % Vandermonde matrix
        Ainv = inv(A);
    else           % tensor product elements
        A = tensorproduct(plocal,porder,pp,dpp);             % Vandermonde matrix
        Ainv = inv(A);
    end 
end

if elemtype==0 % simplex elements
    [nf,nfx,nfy,nfz]=koornwinder(pts,porder,pp,dpp);   % orthogonal shape functions
else           % tensor product elements
    [nf,nfx,nfy,nfz]=tensorproduct(pts,porder,pp,dpp); % orthogonal shape functions
end

% Divide orthogonal shape functions by the Vandemonde matrix to
% obtain nodal shape functions
dim = size(plocal,2);
switch dim
    case 1
        nfs=[nf;nfx] * Ainv;
        nfs=reshape(nfs,[size(nf,1),2,size(nf,2)]);
        nfs=permute(nfs,[3,1,2]);
    case 2
        nfs=[nf;nfx;nfy] * Ainv;
        nfs=reshape(nfs,[size(nf,1),3,size(nf,2)]);
        nfs=permute(nfs,[3,1,2]);
    case 3
        nfs=[nf;nfx;nfy;nfz] * Ainv;
        nfs=reshape(nfs,[size(nf,1),4,size(nf,2)]);
        nfs=permute(nfs,[3,1,2]);
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end

end



function [f,fx,fy,fz] = tensorproduct(x,porder,pp,dpp)
%TENSORPRODUCT Vandermonde matrices for tensor product polynomials in [0,1]^d
%
%   [F,FX,FY,FZ] = TENSORPRODUCT(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+1)*(PORDER+1)
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%      FZ:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. z (npoints,npoly)
%

flag = 0;
porderScalar = porder;

[n,dim] = size(x);

switch (dim)
    case 1 % 1D 
        [f,fx] = koornwinder(x,porder,pp,dpp);  % Legendre basis
        fy = [];
        fz = [];
    case 2 % 2D
        if length(porder)==1
            porder = [porder porder];
        end
        [g1,gx] = koornwinder(x(:,1),porder(1),pp,dpp); % Legendre basis in x direction     
        [g2,gy] = koornwinder(x(:,2),porder(2),pp,dpp); % Legendre basis in y direction         
        f  = zeros(n,prod(porder+1));
        fx = f;
        fy = f;        
        fz = [];
        % perform tensor product to obtain the shape functions and their 
        % derivatives on the unit square
        for ii=1:n
            f(ii,:) =  kron(g2(ii,:),g1(ii,:));
            fx(ii,:) = kron(g2(ii,:),gx(ii,:));
            fy(ii,:) = kron(gy(ii,:),g1(ii,:));
        end
    case 3
        if length(porder)==1
            porder=[porder porder porder];
        end
        
        if flag == 0
            [g1,gx]=koornwinder(x(:,1),porder(1),pp,dpp); % Legendre basis in x direction         
            [g2,gy]=koornwinder(x(:,2),porder(2),pp,dpp); % Legendre basis in y direction             
            [g3,gz]=koornwinder(x(:,3),porder(3),pp,dpp); % Legendre basis in z direction  
        elseif flag == 1
            var1 = repmat(0:porderScalar,[n,1]);
            var2 = zeros(n,1);
            var3 = repmat((0:porderScalar-1),[n,1]);
            aux = [repmat(x(:,1),[1,porderScalar+1]).^var1, var2, (x(:,1)*(1:porderScalar)).^var3];
            g1 = aux(:,1:porderScalar+1);
            gx = aux(:,porderScalar+2:end);
            aux = [repmat(x(:,2),[1,porderScalar+1]).^var1, var2, (x(:,2)*(1:porderScalar)).^var3];
            g2 = aux(:,1:porderScalar+1);
            gy = aux(:,porderScalar+2:end);
            aux = [repmat(x(:,3),[1,porderScalar+1]).^var1, var2, (x(:,3)*(1:porderScalar)).^var3];
            g3 = aux(:,1:porderScalar+1);
            gz = aux(:,porderScalar+2:end);
        end
        
        f  = zeros(n,prod(porder+1));
        fx = f;
        fy = f;
        fz = f;
        % perform tensor product to obtain the shape functions and their 
        % derivatives on the unit cube
        for ii=1:n
            f(ii,:) =  kron(g3(ii,:),kron(g2(ii,:),g1(ii,:)));
            fx(ii,:) = kron(g3(ii,:),kron(g2(ii,:),gx(ii,:)));
            fy(ii,:) = kron(g3(ii,:),kron(gy(ii,:),g1(ii,:)));
            fz(ii,:) = kron(gz(ii,:),kron(g2(ii,:),g1(ii,:)));
        end        
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end

end


function [f,fx,fy,fz] = koornwinder(x,porder,pp,dpp)
%KOORNWINDER Vandermonde matrices for Koornwinder polynomials on simplex elements. 
%            In 1D the Legendre polynomilas normalized to [0,1] are used. 
%            In 2D the koornwinder basis is normalized to the (0,0),(1,0),(0,1) triangle. 
%            In 3D the koornwinder basis is normalized to the 
%            (0,0,0),(1,0,0),(0,1,0),(0,0,1) tetrahedron. The basis are orthonormal.
%
%   [F,FX,FY,FZ]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+2)/2
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%      FZ:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. z (npoints,npoly)
%

dim=size(x,2);

switch dim
 case 1
     [f,fx]=koornwinder1d(x,porder,pp,dpp);      % 1D 
     fy=[];
     fz=[];
 case 2
     [f,fx,fy]=koornwinder2d(x,porder);   % 2D
     fz=[];
 case 3
     [f,fx,fy,fz]=koornwinder3d(x,porder);% 3D           
 otherwise
     error('Only can handle dim=1, dim=2 or dim=3');
end

end


function [f,fx] = koornwinder1d(x,p,pp,dpp)
%KOORNWINDER1D Vandermonde matrix for Legenedre polynomials in [0,1]
%   [F,FX]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (NPOINTS)
%      P:         Maximum order of the polynomials consider. That
%                 is all polynomials of degree up to P, NPOLY=P+1
%      F:         Vandermonde matrix (NPOINTS,NPOLY)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (NPOINTS,NPOLY)
%

x = 2*x-1;

f=zeros(numel(x),p+1);
fx=zeros(numel(x),p+1);

for ii=0:p
    pval  = polyval(pp{ii+1},x);
    dpval = polyval(dpp{ii+1},x);
    f(:,ii+1) = pval;
    fx(:,ii+1) = dpval;
end

fx = 2*fx;

end


function [f,fx,fy] = koornwinder2d(x,p)
%KOORNWINDER2D Vandermonde matrix for Koornwinder polynomials in 
%              the master triangle [0,0]-[1,0]-[0,1]
%   [F,FX,FY]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+2)/2
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%

x = 2*x-1;
npol=prod((double(p)+(1:2))./(1:2));
f=zeros(size(x,1),npol);
fx=zeros(size(x,1),npol);
fy=zeros(size(x,1),npol);

pq = pascalindex2d(npol);

xc = x;
xc(:,2) = min( 0.99999999, xc(:,2)); % avoid singularity

e(:,1) = 2*(1+xc(:,1))./(1-xc(:,2))-1;
e(:,2) = xc(:,2);
ii = find(x(:,2) == 1);
% Correct values for function evaluation
e(ii,1) = -1;
e(ii,2) =  1;

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end

    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));
    
    % Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1));

    f(:,ii) = fc*pval.*qval;
end


% Use displaced coordinate for derivative evaluation
e(:,1) = 2*(1+xc(:,1))./(1-xc(:,2))-1;
e(:,2) = xc(:,2);
de1(:,1) = 2./(1-xc(:,2));
de1(:,2) = 2*(1+xc(:,1))./(1-xc(:,2)).^2;

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end

    dpp = polyder(pp);
    dqp = polyder(qp);

    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));

    dpval = polyval(dpp,e(:,1));
    dqval = polyval(dqp,e(:,2));
    
    % Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1));

    fx(:,ii) = fc*dpval.*qval.*de1(:,1);
    fy(:,ii) = fc*(dpval.*qval.*de1(:,2) + pval.*dqval);
end
fx = 2*fx;
fy = 2*fy;

end


function pq = pascalindex2d(p)
l=1;
for i=0:p
    for j=0:i
        pq(l,1)=i-j;
        pq(l,2)=j;
        l = l+1;
        if l>p 
           return;
        end
    end
end

end


function [f,fx,fy,fz] = koornwinder3d(x,p)
%KOORNWINDER2D Vandermonde matrix for Koornwinder polynomials in 
%              the master tetrahedron [0,0,0]-[1,0,0]-[0,1,0]-[0,0,1]
%   [F,FX,FY]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points where the polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+2)*(PORDER+3)/6
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%      FZ:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. Z (npoints,npoly)
%

x = 2*x-1;
npol=(p+1)*(p+2)*(p+3)/6;
f=zeros(size(x,1),npol);
fx=zeros(size(x,1),npol);
fy=zeros(size(x,1),npol);
fz=zeros(size(x,1),npol);

pq = pascalindex3d(npol);

e  = x;
for ii=1:size(x,1)
    if (x(ii,2)+x(ii,3) ~= 0)
        e(ii,1) = -2*(1+x(ii,1))./(x(ii,2)+x(ii,3))-1;
    else
        e(ii,1) = -1;
    end
    if x(ii,3) ~= 1
        e(ii,2) = 2*(1+x(ii,2))./(1-x(ii,3))-1;
    else
        e(ii,2) = -1;
    end        
    e(ii,3) = x(ii,3);    
end

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    rp = jacobi(pq(ii,3),2*pq(ii,1)+2*pq(ii,2)+2,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end
    for i=1:pq(ii,1)+pq(ii,2)        
        rp = conv([-0.5,0.5],rp);        
    end
    
    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));
    rval = polyval(rp,e(:,3));
    
    % Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1)*2*(pq(ii,1)+pq(ii,2)+pq(ii,3)+2));

    f(:,ii) = fc*pval.*qval.*rval;
end

% Use displaced coordinate for derivative evaluation
xc = x;
for ii=1:size(x,1)
    if (x(ii,2)+x(ii,3) == 0)        
        xc(ii,3) = -1e-8-xc(ii,2);
    end
    if x(ii,3) == 1
        xc(ii,3) = 0.99999999;
    end                
end

e(:,1) = -2*(1+xc(:,1))./(xc(:,2)+xc(:,3))-1;
e(:,2) = 2*(1+xc(:,2))./(1-xc(:,3))-1;
e(:,3) = xc(:,3);
de1(:,1) = -2./(xc(:,2)+xc(:,3));
de1(:,2) = 2*(1+xc(:,1))./(xc(:,2)+xc(:,3)).^2;
de1(:,3) = 2*(1+xc(:,1))./(xc(:,2)+xc(:,3)).^2;
de2(:,1) = 0*xc(:,1); 
de2(:,2) = 2./(1-xc(:,3));
de2(:,3) = 2*(1+xc(:,2))./(1-xc(:,3)).^2;

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    rp = jacobi(pq(ii,3),2*pq(ii,1)+2*pq(ii,2)+2,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end
    for i=1:pq(ii,1)+pq(ii,2)        
        rp = conv([-0.5,0.5],rp);        
    end

    dpp = polyder(pp);
    dqp = polyder(qp);
    drp = polyder(rp);

    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));
    rval = polyval(rp,e(:,3));

    dpval = polyval(dpp,e(:,1));
    dqval = polyval(dqp,e(:,2));
    drval = polyval(drp,e(:,3));
    
    % Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1)*2*(pq(ii,1)+pq(ii,2)+pq(ii,3)+2));

    fx(:,ii) = fc*(dpval.*qval.*rval.*de1(:,1));
    fy(:,ii) = fc*(dpval.*qval.*rval.*de1(:,2) + pval.*dqval.*rval.*de2(:,2));
    fz(:,ii) = fc*(dpval.*qval.*rval.*de1(:,3) + pval.*dqval.*rval.*de2(:,3) + pval.*qval.*drval);
end
fx = 2*fx;
fy = 2*fy;
fz = 2*fz;

end


function pq = pascalindex3d(p)

l=1;
for i=0:p
    for j=0:i
        for k=0:j
            pq(l,1)=i-j;
            pq(l,2)=j-k;
            pq(l,3)=k;
            l = l+1;
            if l>p 
               return;
            end
        end
    end
end

end
