function submesh_index = find_submesh_index(mesh,meshAll)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND_SUBMESH_INDEX Relates the sub implicit or explicit 
% mesh to the full domain mesh
%
%     SUBMESH_INDEX = FIND_SUBMESH_INDEX(MESH, MESHALL)
%
%         MESH      Submesh data structure
%         MESHALL   Full mesh data structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Avg x and y values for each element 
dg_avg_sub = mean(mesh.dgnodes,1);
dg_avg = mean(meshAll.dgnodes,1);

% Number of significant figures to keep
tol = 6;

% Round x and y values
dg_avg = round(dg_avg,tol);
dg_avg_sub = round(dg_avg_sub,tol);

% Simplify to rows of x y
dg_avg = squeeze(dg_avg)';
dg_avg_sub = squeeze(dg_avg_sub)';

% Find corresponding full mesh indices
[~,~,ib] = intersect(dg_avg_sub, dg_avg,'rows','stable');

% Full mesh indices corresponding to submesh
submesh_index = ib;

end