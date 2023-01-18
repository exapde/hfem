
function mesh = get_hField(mesh, limitDGheights)

% This function calculates the mesh size associated with each DG node.
% Ref. David Moro PhD thesis.
%
% Observations:
% - Works properly for non-curved elements. Not confirmed with curved 
%   elements (another approach could be employed then).

if nargin < 2; limitDGheights = 0; end

numDim = mesh.nd;

if numDim == 2
    mesh = get_hField2d(mesh, limitDGheights);
elseif numDim == 3
    mesh = get_hField3d(mesh, limitDGheights);
else
	error('Number of dimensions not implemented.')
end
    
end



function mesh = get_hField2d(mesh, limitDGheights)

% Author: Pablo Fernandez
% Date: Sat, Sep 13, 2014

if nargin < 2; limitDGheights = 0; end

npv = size(mesh.dgnodes,1);

% Creation of a matrix that contains the node to element connectivity:
number_elements_adjacent_to_node = zeros(mesh.np,1);
node_to_elem = zeros(max(mesh.t(:)),20);      % We assume at least 20 elements are connected to a node. It is ok if this number is too large or too small.
for cont_elem=1:mesh.ne
    nodes_in_elem = mesh.t(cont_elem,:);
    for nodes_in_elem_runner = nodes_in_elem
        number_elements_adjacent_to_node(nodes_in_elem_runner) = number_elements_adjacent_to_node(nodes_in_elem_runner)+1;
        node_to_elem(nodes_in_elem_runner,number_elements_adjacent_to_node(nodes_in_elem_runner)) = cont_elem;
    end
end
node_to_elem = node_to_elem(:,1:max(number_elements_adjacent_to_node(:)));


% STEP #1:
% Calculation of the characteristic size of each element (minimum of the
% element heights):
h_elem = zeros(mesh.ne,1);
if mesh.elemtype == 0
    for cont_elem=1:mesh.ne

        nodes_in_elem = mesh.t(cont_elem,:);

        coor_node_1 = mesh.p(nodes_in_elem(1),:);
        coor_node_2 = mesh.p(nodes_in_elem(2),:);
        coor_node_3 = mesh.p(nodes_in_elem(3),:);

        aux_n_face_12 = (coor_node_2-coor_node_1)/norm((coor_node_2-coor_node_1));
        n_face_12(1) = -aux_n_face_12(2);
        n_face_12(2) = aux_n_face_12(1);

        aux_n_face_23 = (coor_node_3-coor_node_2)/norm((coor_node_3-coor_node_2));
        n_face_23(1) = -aux_n_face_23(2);
        n_face_23(2) = aux_n_face_23(1);

        aux_n_face_31 = (coor_node_1-coor_node_3)/norm((coor_node_1-coor_node_3));
        n_face_31(1) = -aux_n_face_31(2);
        n_face_31(2) = aux_n_face_31(1);

        h_elem(cont_elem) = min([abs(n_face_23*(coor_node_2-coor_node_1)'),abs(n_face_31*(coor_node_3-coor_node_2)'),abs(n_face_12*(coor_node_1-coor_node_3)')]);
    end
elseif mesh.elemtype == 1
    for cont_elem=1:mesh.ne

        nodes_in_elem = mesh.t(cont_elem,:);

        coor_node_1 = mesh.p(nodes_in_elem(1),:);
        coor_node_2 = mesh.p(nodes_in_elem(2),:);
        coor_node_3 = mesh.p(nodes_in_elem(3),:);
        coor_node_4 = mesh.p(nodes_in_elem(4),:);

        h_elem(cont_elem) = min([norm(coor_node_2-coor_node_1),norm(coor_node_3-coor_node_1),norm(coor_node_4-coor_node_3),norm(coor_node_1-coor_node_4),norm(coor_node_3-coor_node_2),norm(coor_node_4-coor_node_2)]);
    end
end


% STEP #2:
% Calculation of the characteristic size of each vertex (average of the
% characteristic size the adjacent elements):
h_node = zeros(mesh.np,1);
for cont_node=1:mesh.np
    if number_elements_adjacent_to_node(cont_node) > 0
        adjacent_elements_to_node = node_to_elem(cont_node,:);
        adjacent_elements_to_node = adjacent_elements_to_node(1:number_elements_adjacent_to_node(cont_node));

        for adjacent_elements_to_node_runner = adjacent_elements_to_node
            h_node(cont_node) = h_node(cont_node) + h_elem(adjacent_elements_to_node_runner);
        end
        h_node(cont_node) = h_node(cont_node)/number_elements_adjacent_to_node(cont_node);
    else
        h_node(cont_node) = 0;
    end
end


% STEP #3:
% Calculation of the characteristic size of each DG node (linear
% inteprolation of the characteristic size of each vertex):
mesh.hField = zeros(length(mesh.dgnodes(:,1,1)),mesh.ne);
xi1 = mesh.plocal(:,1);
xi2 = mesh.plocal(:,2);
if mesh.elemtype == 0
    for cont_e = 1:mesh.ne
%         p1 = mesh.p(mesh.t(cont_e,1),:);
%         p2 = mesh.p(mesh.t(cont_e,2),:);
%         p3 = mesh.p(mesh.t(cont_e,3),:);

        h1 = h_node(mesh.t(cont_e,1));
        h2 = h_node(mesh.t(cont_e,2));
        h3 = h_node(mesh.t(cont_e,3));
        hMin = min([h1,h2,h3]);
        hMax = max([h1,h2,h3]);

%         A = [h1 h2 h3]/[[p1';1] [p2';1] [p3';1]]; % Matrix with coefficient for linear interpolation.
% 
%         for cont_DG_node_inside_elem = 1:length(mesh.dgnodes(:,1,cont_elem))
%             aux_hField = A*[(mesh.dgnodes(cont_DG_node_inside_elem,1:2,cont_e))';1];
% 
%             % Approach to mitigate problems with pathologically curved
%             % elements:
%             if limitDGheights == 1
%                 aux_hField = max([aux_hField,hMin]);
%                 aux_hField = min([aux_hField,hMax]);
%             end
% 
%             mesh.hField(cont_DG_node_inside_elem,cont_e) = aux_hField;
%         end
        
        A = [h1 h2 h3]/[[0;0;1] [1;0;1] [0;1;1]]; % Matrix with coefficient for linear interpolation.
        aux_hField = A*[xi1(:)'; xi2(:)'; ones(1,npv)];

        % Approach to mitigate problems with pathologically curved
        % elements:
        if limitDGheights == 1
            aux_hField = max([aux_hField,repmat(hMin,[npv,1])],[],2);
            aux_hField = min([aux_hField,repmat(hMax,[npv,1])],[],2);
        end
        
        mesh.hField(:,cont_e) = aux_hField;
    end
elseif mesh.elemtype == 1
    for cont_e = 1:mesh.ne
%         p1 = mesh.p(mesh.t(cont_e,1),:);
%         p2 = mesh.p(mesh.t(cont_e,2),:);
%         p3 = mesh.p(mesh.t(cont_e,3),:);
%         p4 = mesh.p(mesh.t(cont_e,4),:);
%         
%         v12 = p2-p1;
%         v43 = p3-p4;
%         v14 = p4-p1;
%         v1p = mesh.dgnodes(:,1:2,cont_e) - repmat(p1(:)',[npv,1]);
%         v4p = mesh.dgnodes(:,1:2,cont_e) - repmat(p4(:)',[npv,1]);
        
%         a12 = v1p*v12(:)/norm(v12)^2;
%         a43 = v4p*v43(:)/norm(v43)^2;
%         a14 = v1p*v14(:)/norm(v14)^2;

        h1 = h_node(mesh.t(cont_e,1));
        h2 = h_node(mesh.t(cont_e,2));
        h3 = h_node(mesh.t(cont_e,3));
        h4 = h_node(mesh.t(cont_e,4));
        hMin = min([h1,h2,h3,h4]);
        hMax = max([h1,h2,h3,h4]);
        
%         h12 = h2*a12 + h1*(1-a12);
%         h43 = h3*a43 + h4*(1-a43);
%         aux_hField = h43.*a14 + h12.*(1-a14);

        h12 = (1-xi1)*h1 + xi1*h2;
        h43 = (1-xi1)*h4 + xi1*h3;
        aux_hField = (1-xi2).*h12 + xi2.*h43;
        
        if limitDGheights == 1
            aux_hField = max([aux_hField,repmat(hMin,[npv,1])],[],2);
            aux_hField = min([aux_hField,repmat(hMax,[npv,1])],[],2);
        end
        
        mesh.hField(:,cont_e) = aux_hField;
    end
end

end


function mesh = get_hField3d(mesh, limitDGheights)

% Author: Pablo Fernandez
% Date: Nov, 04, 2014 (adapted from 2D version)

if nargin < 2; limitDGheights = 0; end

npv = size(mesh.dgnodes,1);

% Creation of a matrix that contains the node to element connectivity:
number_elements_adjacent_to_node = zeros(mesh.np,1);
node_to_elem = zeros(max(mesh.t(:)),20);      % We assume at least 20 elements are connected to a node. It is ok if this number is too large or too small.
for cont_elem=1:mesh.ne
    nodes_in_elem = mesh.t(cont_elem,:);
    for nodes_in_elem_runner = nodes_in_elem
        number_elements_adjacent_to_node(nodes_in_elem_runner) = number_elements_adjacent_to_node(nodes_in_elem_runner)+1;
        node_to_elem(nodes_in_elem_runner,number_elements_adjacent_to_node(nodes_in_elem_runner)) = cont_elem;
    end
end
node_to_elem = node_to_elem(:,1:max(number_elements_adjacent_to_node(:)));


% STEP #1:
% Calculation of the characteristic size of each element (minimum of the
% element heights):
h_elem = zeros(mesh.ne,1);
if mesh.elemtype == 0
    for cont_elem=1:mesh.ne

        nodes_in_elem = mesh.t(cont_elem,:);

        coor_node_1 = mesh.p(nodes_in_elem(1),:);
        coor_node_2 = mesh.p(nodes_in_elem(2),:);
        coor_node_3 = mesh.p(nodes_in_elem(3),:);
        coor_node_4 = mesh.p(nodes_in_elem(4),:);

        vector12 = (coor_node_2-coor_node_1)/norm((coor_node_2-coor_node_1));
        vector13 = (coor_node_3-coor_node_1)/norm((coor_node_3-coor_node_1));
        vector14 = (coor_node_4-coor_node_1)/norm((coor_node_4-coor_node_1));
        vector23 = (coor_node_3-coor_node_2)/norm((coor_node_3-coor_node_2));
        vector24 = (coor_node_4-coor_node_2)/norm((coor_node_4-coor_node_2));
        % vector34 = (coor_node_4-coor_node_3)/norm((coor_node_4-coor_node_3));    <- Not necessary.

        n_face_123 = cross(vector12,vector13);
        n_face_123 = n_face_123/norm(n_face_123);

        n_face_134 = cross(vector13,vector14);
        n_face_134 = n_face_134/norm(n_face_134);

        n_face_124 = cross(vector12,vector14);
        n_face_124 = n_face_124/norm(n_face_124);

        n_face_234 = cross(vector23,vector24);
        n_face_234 = n_face_234/norm(n_face_234);

        h_elem(cont_elem) = min([abs(n_face_123*(coor_node_4-coor_node_1)'), ...
            abs(n_face_134*(coor_node_2-coor_node_1)'), ...
            abs(n_face_124*(coor_node_3-coor_node_1)'), ...
            abs(n_face_234*(coor_node_1-coor_node_2)')]);
    end
elseif mesh.elemtype == 1
    for cont_elem=1:mesh.ne

        nodes_in_elem = mesh.t(cont_elem,:);

        coor_node_1 = mesh.p(nodes_in_elem(1),:);
        coor_node_2 = mesh.p(nodes_in_elem(2),:);
        coor_node_3 = mesh.p(nodes_in_elem(3),:);
        coor_node_4 = mesh.p(nodes_in_elem(4),:);
        coor_node_5 = mesh.p(nodes_in_elem(5),:);
        coor_node_6 = mesh.p(nodes_in_elem(6),:);
        coor_node_7 = mesh.p(nodes_in_elem(7),:);
        coor_node_8 = mesh.p(nodes_in_elem(8),:);

        vector12 = norm((coor_node_2-coor_node_1));
        vector13 = norm((coor_node_3-coor_node_1));
        vector14 = norm((coor_node_4-coor_node_1));
        vector15 = norm((coor_node_5-coor_node_1));
        vector16 = norm((coor_node_6-coor_node_1));
        vector17 = norm((coor_node_7-coor_node_1));
        vector18 = norm((coor_node_8-coor_node_1));
        vector23 = norm((coor_node_3-coor_node_2));
        vector24 = norm((coor_node_4-coor_node_2));
        vector25 = norm((coor_node_5-coor_node_2));
        vector26 = norm((coor_node_6-coor_node_2));
        vector27 = norm((coor_node_7-coor_node_2));
        vector28 = norm((coor_node_8-coor_node_2));
        vector34 = norm((coor_node_4-coor_node_3));
        vector35 = norm((coor_node_5-coor_node_3));
        vector36 = norm((coor_node_6-coor_node_3));
        vector37 = norm((coor_node_7-coor_node_3));
        vector38 = norm((coor_node_8-coor_node_3));
        vector45 = norm((coor_node_5-coor_node_4));
        vector46 = norm((coor_node_6-coor_node_4));
        vector47 = norm((coor_node_7-coor_node_4));
        vector48 = norm((coor_node_8-coor_node_4));
        vector56 = norm((coor_node_6-coor_node_5));
        vector57 = norm((coor_node_7-coor_node_5));
        vector58 = norm((coor_node_8-coor_node_5));
        vector67 = norm((coor_node_7-coor_node_6));
        vector68 = norm((coor_node_8-coor_node_6));
        vector78 = norm((coor_node_8-coor_node_7));

        h_elem(cont_elem) = min([vector12,vector13,vector14,vector15,vector16,...
            vector17,vector18,vector23,vector24,vector25,vector26,vector27,...
            vector28,vector34,vector35,vector36,vector37,vector38,vector45,...
            vector46,vector47,vector48,vector56,vector57,vector58,vector67,...
            vector68,vector78]);
    end
end


% STEP #2:
% Calculation of the characteristic size of each vertex (average of the
% characteristic size the adjacent elements):
h_node = zeros(mesh.np,1);
for cont_node=1:mesh.np
    if number_elements_adjacent_to_node(cont_node) > 0
        adjacent_elements_to_node = node_to_elem(cont_node,:);
        adjacent_elements_to_node = adjacent_elements_to_node(1:number_elements_adjacent_to_node(cont_node));

        for adjacent_elements_to_node_runner = adjacent_elements_to_node
            h_node(cont_node) = h_node(cont_node) + h_elem(adjacent_elements_to_node_runner);
        end
        h_node(cont_node) = h_node(cont_node)/number_elements_adjacent_to_node(cont_node);
    else
        h_node(cont_node) = 0;
    end
end


% STEP #3:
% Calculation of the characteristic size of each DG node (linear
% inteprolation of the characteristic size of each vertex):
mesh.hField = zeros(length(mesh.dgnodes(:,1,1)),mesh.ne);
xi1 = mesh.plocal(:,1);
xi2 = mesh.plocal(:,2);
xi3 = mesh.plocal(:,3);
if mesh.elemtype == 0
    for cont_e = 1:mesh.ne
%         p1 = mesh.p(mesh.t(cont_e,1),:);
%         p2 = mesh.p(mesh.t(cont_e,2),:);
%         p3 = mesh.p(mesh.t(cont_e,3),:);
%         p4 = mesh.p(mesh.t(cont_e,4),:);

        h1 = h_node(mesh.t(cont_e,1));
        h2 = h_node(mesh.t(cont_e,2));
        h3 = h_node(mesh.t(cont_e,3));
        h4 = h_node(mesh.t(cont_e,4));
        hMin = min([h1,h2,h3,h4]);
        hMax = max([h1,h2,h3,h4]);

%         A = [h1 h2 h3 h4]/[[p1';1] [p2';1] [p3';1] [p4';1]]; % Matrix with coefficient for linear interpolation.
% 
%         for cont_DG_node_inside_elem = 1:length(mesh.dgnodes(:,1,cont_elem))
%             aux_hField = A*[(mesh.dgnodes(cont_DG_node_inside_elem,1:3,cont_e))';1];
% 
%             % Approach to mitigate problems with pathologically curved
%             % elements:
%             if limitDGheights == 1
%                 aux_hField = max([aux_hField,hMin]);
%                 aux_hField = min([aux_hField,hMax]);
%             end
%             
%             mesh.hField(cont_DG_node_inside_elem,cont_e) = aux_hField;
%         end
        
        A = [h1 h2 h3 h4]/[[0;0;0;1] [1;0;0;1] [0;1;0;1] [0;0;1;1]]; % Matrix with coefficient for linear interpolation.
        aux_hField = A*[xi1(:)'; xi2(:)'; xi3(:)'; ones(1,npv)];

        % Approach to mitigate problems with pathologically curved
        % elements:
        if limitDGheights == 1
            aux_hField = max([aux_hField,repmat(hMin,[npv,1])],[],2);
            aux_hField = min([aux_hField,repmat(hMax,[npv,1])],[],2);
        end
        
        mesh.hField(:,cont_e) = aux_hField;
    end
elseif mesh.elemtype == 1
    for cont_e = 1:mesh.ne
%         p1 = mesh.p(mesh.t(cont_e,1),:);
%         p2 = mesh.p(mesh.t(cont_e,2),:);
%         p3 = mesh.p(mesh.t(cont_e,3),:);
%         p4 = mesh.p(mesh.t(cont_e,4),:);
%         p5 = mesh.p(mesh.t(cont_e,5),:);
%         p6 = mesh.p(mesh.t(cont_e,6),:);
%         p7 = mesh.p(mesh.t(cont_e,7),:);
%         p8 = mesh.p(mesh.t(cont_e,8),:);
        
%         v12 = p2-p1;
%         v43 = p3-p4;
%         v14 = p4-p1;
%         v56 = p6-p5;
%         v87 = p7-p8;
%         v58 = p8-p5;
%         v15 = p5-p1;
%         v1p = mesh.dgnodes(:,1:3,cont_e) - repmat(p1(:)',[npv,1]);
%         v4p = mesh.dgnodes(:,1:3,cont_e) - repmat(p4(:)',[npv,1]);
%         v5p = mesh.dgnodes(:,1:3,cont_e) - repmat(p5(:)',[npv,1]);
%         v8p = mesh.dgnodes(:,1:3,cont_e) - repmat(p8(:)',[npv,1]);
        
%         a12 = v1p*v12(:)/norm(v12)^2;
%         a43 = v4p*v43(:)/norm(v43)^2;
%         a56 = v5p*v56(:)/norm(v56)^2;
%         a87 = v8p*v87(:)/norm(v87)^2;
%         a14 = v1p*v14(:)/norm(v14)^2;
%         a58 = v5p*v58(:)/norm(v58)^2;
%         a15 = v1p*v15(:)/norm(v15)^2;

        h1 = h_node(mesh.t(cont_e,1));
        h2 = h_node(mesh.t(cont_e,2));
        h3 = h_node(mesh.t(cont_e,3));
        h4 = h_node(mesh.t(cont_e,4));
        h5 = h_node(mesh.t(cont_e,5));
        h6 = h_node(mesh.t(cont_e,6));
        h7 = h_node(mesh.t(cont_e,7));
        h8 = h_node(mesh.t(cont_e,8));
        hMin = min([h1,h2,h3,h4,h5,h6,h7,h8]);
        hMax = max([h1,h2,h3,h4,h5,h6,h7,h8]);
        
%         h12 = h2*a12 + h1*(1-a12);
%         h43 = h3*a43 + h4*(1-a43);
%         h56 = h6*a56 + h5*(1-a56);
%         h87 = h7*a87 + h8*(1-a87);
%         h1234 = h43.*a14 + h12.*(1-a14);
%         h5678 = h87.*a58 + h56.*(1-a58);
%         aux_hField = h5678.*a15 + h1234.*(1-a15);
        
        h12 = (1-xi1)*h1 + xi1*h2;
        h43 = (1-xi1)*h4 + xi1*h3;
        h1234 = (1-xi2).*h12 + xi2.*h43;

        h56 = (1-xi1)*h5 + xi1*h6;
        h87 = (1-xi1)*h8 + xi1*h7;
        h5678 = (1-xi2).*h56 + xi2.*h87;

        aux_hField = (1-xi3).*h1234 + xi3.*h5678;
        
        if limitDGheights == 1
            aux_hField = max([aux_hField,repmat(hMin,[npv,1])],[],2);
            aux_hField = min([aux_hField,repmat(hMax,[npv,1])],[],2);
        end
        
        mesh.hField(:,cont_e) = aux_hField;
    end 
end

end


% function Px = oblProj(Range,Null,x)
% 
% dim = size(Range,2);
% dimRange = size(Range,1);
% dimNull = size(Null,1);
% 
% if dim ~= size(Null,2); error('Dimensions do not match.'); end
% if dimRange+dimNull ~= dim; error('Dimensions do not match.'); end
% if length(x) ~= dim; error('Dimensions do not match.'); end
% 
% eigMatrix = [Range', Null'];
% P = eigMatrix * diag([ones(1,dimRange), zeros(1,dimNull)]) / eigMatrix;
% Px = P*x(:);
% 
% end
