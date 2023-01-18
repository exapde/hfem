
function p1CG = p0DG_to_p1CG(p0DG, mesh, limitp1CGvalues)

% Observations:
% - Works properly for non-curved elements. Not confirmed with curved 
%   elements (another approach could be employed then).

if nargin < 3; limitp1CGvalues = 0; end
if length(p0DG) ~= mesh.ne; error('p0DG input in p0DG_to_p1CG.m must be a vector of length mesh.ne.'); end

numDim = mesh.nd;

if numDim == 2
    p1CG = p0DG_to_p1CG_2d(p0DG, mesh, limitp1CGvalues);
elseif numDim == 3
    p1CG = p0DG_to_p1CG_3d(p0DG, mesh, limitp1CGvalues);
else
	error('Number of dimensions not implemented.')
end
    
end


function p1CG = p0DG_to_p1CG_2d(p0DG, mesh, limitp1CGvalues)

if nargin < 3; limitp1CGvalues = 0; end

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
% Calculation of the characteristic size of each vertex (average of the
% characteristic size the adjacent elements):
field_node = zeros(mesh.np,1);
for cont_node=1:mesh.np
    if number_elements_adjacent_to_node(cont_node) > 0
        adjacent_elements_to_node = node_to_elem(cont_node,:);
        adjacent_elements_to_node = adjacent_elements_to_node(1:number_elements_adjacent_to_node(cont_node));

        for adjacent_elements_to_node_runner = adjacent_elements_to_node
            field_node(cont_node) = field_node(cont_node) + p0DG(adjacent_elements_to_node_runner);
        end
        field_node(cont_node) = field_node(cont_node)/number_elements_adjacent_to_node(cont_node);
    else
        field_node(cont_node) = 0;
    end
end


% STEP #2:
% Calculation of the characteristic size of each DG node (linear
% inteprolation of the characteristic size of each vertex):
p1CG = zeros(length(mesh.dgnodes(:,1,1)),mesh.ne);
xi1 = mesh.plocal(:,1);
xi2 = mesh.plocal(:,2);
if mesh.elemtype == 0
    for cont_e = 1:mesh.ne
%         p1 = mesh.p(mesh.t(cont_e,1),:);
%         p2 = mesh.p(mesh.t(cont_e,2),:);
%         p3 = mesh.p(mesh.t(cont_e,3),:);

        field1 = field_node(mesh.t(cont_e,1));
        field2 = field_node(mesh.t(cont_e,2));
        field3 = field_node(mesh.t(cont_e,3));
        fieldMin = min([field1,field2,field3]);
        fieldMax = max([field1,field2,field3]);

%         A = [h1 h2 h3]/[[p1';1] [p2';1] [p3';1]]; % Matrix with coefficient for linear interpolation.
% 
%         for cont_DG_node_inside_elem = 1:length(mesh.dgnodes(:,1,cont_elem))
%             aux_p1CG = A*[(mesh.dgnodes(cont_DG_node_inside_elem,1:2,cont_e))';1];
% 
%             % Approach to mitigate problems with pathologically curved
%             % elements:
%             if limitDGheights == 1
%                 aux_p1CG = max([aux_p1CG,hMin]);
%                 aux_p1CG = min([aux_p1CG,hMax]);
%             end
% 
%             p1CG(cont_DG_node_inside_elem,cont_e) = aux_p1CG;
%         end
        
        A = [field1 field2 field3]/[[0;0;1] [1;0;1] [0;1;1]]; % Matrix with coefficient for linear interpolation.
        aux_p1CG = A*[xi1(:)'; xi2(:)'; ones(1,npv)];

        % Approach to mitigate problems with pathologically curved
        % elements:
        if limitp1CGvalues == 1
            aux_p1CG = max([aux_p1CG,repmat(fieldMin,[npv,1])],[],2);
            aux_p1CG = min([aux_p1CG,repmat(fieldMax,[npv,1])],[],2);
        end
        
        p1CG(:,cont_e) = aux_p1CG;
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

        field1 = field_node(mesh.t(cont_e,1));
        field2 = field_node(mesh.t(cont_e,2));
        field3 = field_node(mesh.t(cont_e,3));
        field4 = field_node(mesh.t(cont_e,4));
        fieldMin = min([field1,field2,field3,field4]);
        fieldMax = max([field1,field2,field3,field4]);
        
%         h12 = h2*a12 + h1*(1-a12);
%         h43 = h3*a43 + h4*(1-a43);
%         aux_p1CG = h43.*a14 + h12.*(1-a14);

        field12 = (1-xi1)*field1 + xi1*field2;
        field43 = (1-xi1)*field4 + xi1*field3;
        aux_p1CG = (1-xi2).*field12 + xi2.*field43;
        
        if limitp1CGvalues == 1
            aux_p1CG = max([aux_p1CG,repmat(fieldMin,[npv,1])],[],2);
            aux_p1CG = min([aux_p1CG,repmat(fieldMax,[npv,1])],[],2);
        end
        
        p1CG(:,cont_e) = aux_p1CG;
    end
end

end


function p1CG = p0DG_to_p1CG_3d(p0DG, mesh, limitp1CGvalues)

if nargin < 3; limitp1CGvalues = 0; end

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
% Calculation of the characteristic size of each vertex (average of the
% characteristic size the adjacent elements):
field_node = zeros(mesh.np,1);
for cont_node=1:mesh.np
    if number_elements_adjacent_to_node(cont_node) > 0
        adjacent_elements_to_node = node_to_elem(cont_node,:);
        adjacent_elements_to_node = adjacent_elements_to_node(1:number_elements_adjacent_to_node(cont_node));

        for adjacent_elements_to_node_runner = adjacent_elements_to_node
            field_node(cont_node) = field_node(cont_node) + p0DG(adjacent_elements_to_node_runner);
        end
        field_node(cont_node) = field_node(cont_node)/number_elements_adjacent_to_node(cont_node);
    else
        field_node(cont_node) = 0;
    end
end


% STEP #2:
% Calculation of the characteristic size of each DG node (linear
% inteprolation of the characteristic size of each vertex):
p1CG = zeros(length(mesh.dgnodes(:,1,1)),mesh.ne);
xi1 = mesh.plocal(:,1);
xi2 = mesh.plocal(:,2);
xi3 = mesh.plocal(:,3);
if mesh.elemtype == 0
    for cont_e = 1:mesh.ne
%         p1 = mesh.p(mesh.t(cont_e,1),:);
%         p2 = mesh.p(mesh.t(cont_e,2),:);
%         p3 = mesh.p(mesh.t(cont_e,3),:);
%         p4 = mesh.p(mesh.t(cont_e,4),:);

        field1 = field_node(mesh.t(cont_e,1));
        field2 = field_node(mesh.t(cont_e,2));
        field3 = field_node(mesh.t(cont_e,3));
        field4 = field_node(mesh.t(cont_e,4));
        fieldMin = min([field1,field2,field3,field4]);
        fieldMax = max([field1,field2,field3,field4]);

%         A = [h1 h2 h3 h4]/[[p1';1] [p2';1] [p3';1] [p4';1]]; % Matrix with coefficient for linear interpolation.
% 
%         for cont_DG_node_inside_elem = 1:length(mesh.dgnodes(:,1,cont_elem))
%             aux_p1CG = A*[(mesh.dgnodes(cont_DG_node_inside_elem,1:3,cont_e))';1];
% 
%             % Approach to mitigate problems with pathologically curved
%             % elements:
%             if limitDGheights == 1
%                 aux_p1CG = max([aux_p1CG,hMin]);
%                 aux_p1CG = min([aux_p1CG,hMax]);
%             end
%             
%             p1CG(cont_DG_node_inside_elem,cont_e) = aux_p1CG;
%         end
        
        A = [field1 field2 field3 field4]/[[0;0;0;1] [1;0;0;1] [0;1;0;1] [0;0;1;1]]; % Matrix with coefficient for linear interpolation.
        aux_p1CG = A*[xi1(:)'; xi2(:)'; xi3(:)'; ones(1,npv)];

        % Approach to mitigate problems with pathologically curved
        % elements:
        if limitp1CGvalues == 1
            aux_p1CG = max([aux_p1CG,repmat(fieldMin,[npv,1])],[],2);
            aux_p1CG = min([aux_p1CG,repmat(fieldMax,[npv,1])],[],2);
        end
        
        p1CG(:,cont_e) = aux_p1CG;
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

        field1 = field_node(mesh.t(cont_e,1));
        field2 = field_node(mesh.t(cont_e,2));
        field3 = field_node(mesh.t(cont_e,3));
        field4 = field_node(mesh.t(cont_e,4));
        field5 = field_node(mesh.t(cont_e,5));
        field6 = field_node(mesh.t(cont_e,6));
        field7 = field_node(mesh.t(cont_e,7));
        field8 = field_node(mesh.t(cont_e,8));
        fieldMin = min([field1,field2,field3,field4,field5,field6,field7,field8]);
        fieldMax = max([field1,field2,field3,field4,field5,field6,field7,field8]);
        
%         h12 = h2*a12 + h1*(1-a12);
%         h43 = h3*a43 + h4*(1-a43);
%         h56 = h6*a56 + h5*(1-a56);
%         h87 = h7*a87 + h8*(1-a87);
%         h1234 = h43.*a14 + h12.*(1-a14);
%         h5678 = h87.*a58 + h56.*(1-a58);
%         aux_p1CG = h5678.*a15 + h1234.*(1-a15);
        
        field12 = (1-xi1)*field1 + xi1*field2;
        field43 = (1-xi1)*field4 + xi1*field3;
        field1234 = (1-xi2).*field12 + xi2.*field43;

        field56 = (1-xi1)*field5 + xi1*field6;
        field87 = (1-xi1)*field8 + xi1*field7;
        field5678 = (1-xi2).*field56 + xi2.*field87;

        aux_p1CG = (1-xi3).*field1234 + xi3.*field5678;
        
        if limitp1CGvalues == 1
            aux_p1CG = max([aux_p1CG,repmat(fieldMin,[npv,1])],[],2);
            aux_p1CG = min([aux_p1CG,repmat(fieldMax,[npv,1])],[],2);
        end
        
        p1CG(:,cont_e) = aux_p1CG;
    end 
end

end
