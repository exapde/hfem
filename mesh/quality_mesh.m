function [ mesh_qual ] = quality_mesh( mesh, master)
%QUALITY_MESH Computes the quality of high order tetrahedra meshes
%   The computation of the quality mesure is based of the Xevi Roca's works
%   mainly "Distortion and Quality Measures for validating and generating
%   high order tetrahedral meshes", International Meshing roundtable, 2013
%
%   MESH      :: mesh structure
%   MASTER    :: master structure
%   MESH_QUAL :: quality for each node and each element of the mesh

% Getting Mesh information
ngv = master.ngv;
npv = master.npv;
nd  = master.nd;
mesh_qual = zeros(mesh.ne,1);

% Vertices positions of the ideal tetrahedron in the reference coordinates
V0 = [0. 0. 1.]';
V1 = [1. 0. 0.]';
V2 = [0. 1. 0.]';
V3 = [1. 1. 1.]';
% Jacobian of the corresponding linear transformation between the Master
% element and the Ideal equilateral tetrahedron
JacMI = [V1-V0 V2-V0 V3-V0];
detJacMI = det(JacMI);
% Jacobian of the inverse transformation (Ideal to Master)
JacIM = inv(JacMI);
% sqrt of the Volume of the ideal tatrahedron
VolI = sqrt(detJacMI/6.);

% Derivatives of shape functions
dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);

nundint = 0;
for i=1:mesh.ne
    
    dg = mesh.dgnodes(:,:,i);
    % compute the Jacobian (determinant) at Gauss points
    JacMP = dshapvt*dg(:,1:nd);
    
    % Transposition...
    JacMP = permute(reshape(JacMP,[ngv nd nd]),[1 3 2]);
    JacMP = reshape(JacMP,[ngv*nd nd]);
    % End Transposition
    
    % Getting the Jacobian of the Ideal to Physical transformation
    JacMP = JacMP*JacIM;

    % Get the determinant and Frobenius norm square
    JacMP = reshape(JacMP,[ngv nd nd]);
    frobN    = frobNorm(JacMP);
    detJacMP = detJac(JacMP);
    detJacMP = 1/2 * (abs(detJacMP) + detJacMP);
    % Check if element is degenerated
    if ~all(detJacMP)
        mesh_qual(i) = 1000.0;
        continue
    end
    
    % Computation of the pointwise distortion measure
    distm = frobN(:) ./ (3*abs(detJacMP(:)).^(2/3));
    
    % Integration over the Master element
    DphiE = detJacMI * master.gwvl' * (distm.^2);
    
    % Check under integration
    if DphiE<0
        DphiE = detJacMI * abs(master.gwvl)' * (distm.^2);
        nundint = nundint+1;
    end
    
    % Normalization by the ideal volume
    mesh_qual(i) = sqrt(DphiE)/VolI;

end

if nundint>0
    fprintf('ATTENTION : %d elements are obviously under-integrated \n', nundint);
end


function [jac] = detJac(Jg)

nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;        
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);        
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);                    
    otherwise
        error('Dimension is not implemented');
end

function [frob] = frobNorm(A)
% Returns the square of the Frobenius Norm of a matrix

nd  = size(A,2);
switch nd
    case 1
        frob = A^2;
    case 2
        frob = A(:,1,1).^2 + A(:,1,2).^2 + A(:,2,1).^2 + A(:,2,2).^2;        
    case 3
        frob = A(:,1,1).^2 + A(:,1,2).^2 + A(:,1,3).^2 + ...
               A(:,2,1).^2 + A(:,2,2).^2 + A(:,2,3).^2 + ...
               A(:,3,1).^2 + A(:,3,2).^2 + A(:,3,3).^2;                    
    otherwise
        error('Dimension is not implemented');
end
