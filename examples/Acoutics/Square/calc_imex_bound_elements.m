function mesh = calc_imex_bound_elements(master, mesh, app, SDG, fc, time, UDG)

bound_ele = mesh.bound_ele;

% Find faces and elements for computations
global_faces = find(bound_ele(:,3) ~= 0);  % search bound_ele matrix for boundary elements
bound_ele = bound_ele(global_faces,:);  % gets rid of the zero rows in bound_ele

bound_ele = cat(2, bound_ele,global_faces); % add in the information about the global faces
sorted_faces = sortrows(bound_ele,3); % sort the matrix by the implicit element number

imp_ele = unique(sorted_faces(:,3)); % find the unique boundary elements
num_imp = length(imp_ele); % number of implicit elements

element_data = zeros(3,5,num_imp);  % allocates structure for element data

% loop over implicit elements
for i = 1:num_imp
    j = imp_ele(i);
    ind = find(bound_ele(:,3) == j);
    
    imex_faces = length(ind);
    
    for k = 1:imex_faces
        element_data(k,:,i) = bound_ele(ind(k),:);
    end
end







end