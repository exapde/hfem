function mesh = create_mesh_imface(mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE_MESH_IMFACE Determines if a face is on an implicit element
%
%    MESH = CREATE_MESH_IMFACE(MESH)
%
%       MESH   Mesh data structure
%       IMFACE   1 for implicit face, 0 for not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nf = size(mesh.f,1);
imface = zeros(nf,1);

for i = 1:nf
    % Determine elements and imex
    e1 = mesh.f(i,3);  % element 1
    e2 = mesh.f(i,4);  % element 2
    
    imex1 = mesh.imex(e1);  
    % Account for domain boundary 
    if e2 > 0
        imex2 = mesh.imex(e2);
    else 
        imex2 = 0;
    end
    
    % Create fimex
    if imex1 == 1 || imex2 == 1
        imface(i) = 1;
    end
end

% Add imface to the mesh data structure
mesh.imface = imface;

end