function mesh = create_mesh_imex(mesh,dt,c,Kmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE_MESH_IMEX Determines if an element should be solved implicitly or
% explicitly based on mesh size
%
%    MESH = CREATE_MESH_IMEX(MESH, DT)
%
%       MESH   Mesh data structure
%       DT     Time step
%       IMEX   1 for implicit, 0 for explicit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tri = size(mesh.t,1);
imex = zeros(tri,1);
p = mesh.porder;

% Find the minimum side length of each element
[he, hf] = meshsize(mesh);
minsize = min(he,[],1);

% Determine if element is implicit or explicit
for i = 1:tri
    if minsize(i)*Kmax < c*dt*p^2
        imex(i) = 1;
    end
end

% Add imex to the mesh data structure
mesh.imex = imex;

end