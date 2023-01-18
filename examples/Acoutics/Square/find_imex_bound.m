function mesh = find_imex_bound(mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND_IMEX_BOUND Determines if a face is a boundary between the implicit
% and explicit regions.  Additionally, lists the neighboring explicit and
% implicit elements and the face number of the boundary face.  
%
%    MESH = FIND_IMEX_BOUND(MESH)
%
%       MESH   Mesh data structure
%       FIMEX   1 for boundary, 0 for non-boundary
%       BOUND_ELE   Column 1: explicit element number
%                   Column 2: explicit element imex boundary local face number
%                   Column 3: implicit element number
%                   Column 4: implicit element imex boundary local face number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nf = size(mesh.f,1);
fimex = zeros(nf,1);
bound_ele = zeros(nf,4);

for i = 1:nf
    % Determine elements and imex
    e1 = mesh.f(i,3);  % element 1
    e2 = mesh.f(i,4);  % element 2
    
    % determines if an element is internal to the entire domain
    if e1 > 0 && e2 > 0
        imex1 = mesh.imex(e1);  
        imex2 = mesh.imex(e2);

        % Create fimex
        if imex1 ~= imex2  % 0 for explicit, 1 for implicit
            fimex(i) = 1;
            
            % Finds local face number
            e1_face = find(mesh.t2f(e1,:) == i);
            e2_face = find(mesh.t2f(e2,:) == i);
            
            % Input into bound_ele
            if imex1 == 0
                bound_ele(i,1) = e1;
                bound_ele(i,2) = e1_face;
                bound_ele(i,3) = e2;
                bound_ele(i,4) = e2_face;
            else
                bound_ele(i,1) = e2;
                bound_ele(i,2) = e2_face;
                bound_ele(i,3) = e1;
                bound_ele(i,4) = e1_face;      
            end
        end
    end
end

% Add to mesh data structure
mesh.fimex = fimex;
mesh.bound_ele = bound_ele;

end