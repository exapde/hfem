function mesh = find_imex_bound(mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND_IMEX_BOUND Determines if a face is a boundary between the implicit
% and explicit regions
%
%    MESH = FIND_IMEX_BOUND(MESH)
%
%       MESH   Mesh data structure
%       FIMEX   1 for boundary, 0 for non-boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nf = size(mesh.f,1);
fimex = zeros(nf,1);

for i = 1:nf
    % Determine elements and imex
    e1 = mesh.f(i,3);  % element 1
    e2 = mesh.f(i,4);  % element 2
    
    if e1 > 0 && e2 > 0
        imex1 = mesh.imex(e1);  
        imex2 = mesh.imex(e2);

        % Create fimex
        if imex1 ~= imex2
            fimex(i) = 1;
        end
    end
end

% Add to mesh data structure
mesh.fimex = fimex;

end