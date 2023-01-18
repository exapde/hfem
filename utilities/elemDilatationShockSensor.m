
function S_j = elemDilatationShockSensor(mesh,master,app,UDG,UH,hFlag,smoothSensorFlag,plotFlag)

if nargin < 6; hFlag = 2; end
if nargin < 7; smoothSensorFlag = 0; end
if nargin < 8; plotFlag = 0; end

k_h = 1.5;

nd = mesh.nd;
ne = mesh.ne;
npv = master.npv;
npf = master.npf;
% nfe = mesh.nfe;
ngv = master.ngv;
ngf = master.ngf;
% porder = mesh.porder;
perm = mesh.perm;
elcon = mesh.elcon;
dgnodes = mesh.dgnodes;
porder = mesh.porder;

shapft = squeeze(master.shapft(:,:,1));
% shapfg = squeeze(master.shapfg(:,:,1));
% shapmf = master.shapmf;
dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
gwfc = master.gwfc;

shapvt = squeeze(master.shapvt(:,:,1));
% shapvg = squeeze(master.shapvg(:,:,1));
shapmv = master.shapmv;
gwvl = master.gwvl;

nch = app.nch;
gam = app.arg{1};
gam1 = gam - 1;

r    = squeeze(UDG(:,1,:));
ru   = squeeze(UDG(:,2,:));
rv   = squeeze(UDG(:,3,:));
if nd == 2
    rw   = zeros(size(UDG,1),size(UDG,3));
    rE   = squeeze(UDG(:,4,:));
elseif nd == 3
    rw   = squeeze(UDG(:,4,:));
    rE   = squeeze(UDG(:,5,:));
end
u    = ru./r;
v    = rv./r;
w    = rw./r;
af    = 0.5*(u.*u+v.*v+w.*w);
p   = gam1*(rE -r.*af);
cst = sqrt((2*gam*p./r + gam1*(u.^2 + v.^2 + w.^2))/(gam+1));

numerator = zeros(ne,1);
denominator = zeros(ne,1);

[~, ~, jacE] = volgeom(shapmv,permute(dgnodes(:,1:nd,:),[1,3,2]));
jacE = reshape(jacE, [ngv,ne]);

UH = reshape(UH, nch, []);
for elem = 1:ne
    facesInElem = mesh.t2f(elem,:);
    facesInElem = facesInElem(:)';
    
    pnE = dgnodes(:,1:nd,elem);
    udgE = UDG(:,:,elem);
    udggE = shapvt*udgE;
    
    localFace = 0;
    for face = facesInElem
        localFace = localFace + 1;
        
        pnF = pnE(perm(:,localFace),:);
        pgF = shapft*pnF;

        udggF = shapft*udgE(perm(:,localFace),:);

        uh = UH(:,elcon(:,elem))';
        uhgF = shapft*uh((localFace-1)*npf+1:localFace*npf,:);

        dpgF = dshapft*pnF;
        dpgF = permute(reshape(dpgF,[ngf nd-1 nd]), [1 3 2]);                

        if nd==2
            jacF   = sqrt(dpgF(:,1).^2+dpgF(:,2).^2);
            nlgF   = [dpgF(:,2),-dpgF(:,1)];
            nlgF   = bsxfun(@rdivide, nlgF, jacF);
        elseif nd==3
            nlgF(:,1) = dpgF(:,2,1).*dpgF(:,3,2) - dpgF(:,3,1).*dpgF(:,2,2);
            nlgF(:,2) = dpgF(:,3,1).*dpgF(:,1,2) - dpgF(:,1,1).*dpgF(:,3,2);
            nlgF(:,3) = dpgF(:,1,1).*dpgF(:,2,2) - dpgF(:,2,1).*dpgF(:,1,2);
            jacF = sqrt(nlgF(:,1).^2+nlgF(:,2).^2+nlgF(:,3).^2);
            nlgF   = bsxfun(@rdivide, nlgF, jacF);
        end            

%         fh = fhat(nlgF,pgF,udggF,uhgF,app.arg,0); 
%         numerator(elem) = numerator(elem) -
%         sum(fh(:,1).*jacF(:).*gwfc(:));       % Flux of momentum
        v_n = sum(udggF(:,2:1+nd).*nlgF(:,1:nd)./repmat(udggF(:,1),[1,nd]),2);
        numerator(elem) = numerator(elem) - sum(v_n(:).*jacF(:).*gwfc(:));     % Flux of velocity 
    end
    
%     denominator(elem) = sum(gwvl.*jacE(:,elem).*udggE(:,1));    % mass in element
    denominator(elem) = sum(gwvl.*jacE(:,elem));                % volume in element
end

% Get element size:
if hFlag == 0       % Use constant h = 0.02
%     hElem = 0.02 * ones(ne,1);
    mesh = getElementSize(mesh);
    hField_g = shapvt*squeeze(mesh.hField);

    hElem = zeros(ne,1);
    for elem=1:ne
        hElem(elem) = sum(gwvl .* hField_g(:,elem) .* jacE(:,elem)) / sum(gwvl .*jacE(:,elem));
    end
    hElem = min(0.02,hElem(:));
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

S_j_p0DG = hElem .* numerator ./ denominator;
S_j_p1CG = (k_h/porder) * p0DG_to_p1CG(S_j_p0DG, mesh) ./ cst;

S_j_p0DG = (k_h/porder) * repmat(S_j_p0DG(:)',[npv,1]) ./ cst;

if smoothSensorFlag == 0
    S_j = S_j_p0DG;
elseif smoothSensorFlag == 1
    S_j = S_j_p1CG;
end

if plotFlag; figure(1); scaplot(mesh,S_j); end

end
