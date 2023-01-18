function udg = udgproj(mesh,mastersubgrid,u)

ns = max(mesh.subgrids);
npm = size(mesh.dgnodes,1);
ne  = size(mesh.dgnodes,3);

for j=1:ns
    nc  = size(u{j},2);
    if nc > 0
        break;
    end
end

udg = zeros(npm,nc,ne);
%mesht = mesh;
for j = 1:ns     
    ind = find(mesh.subgrids==j);    
    if isempty(ind)==0
        if j==1
            udg(:,:,ind) = u{j};
        else
            nj = length(ind);
            ng = size(mastersubgrid{j}.elemnodes,3);        
            npv = size(u{j},1);
            %mesht.dgnodes = mesh.dgnodes(:,:,ind);                                             
            udg(:,:,ind) = subgridprojection(mastersubgrid{j},mastersubgrid{1},mesh.dgnodes(:,:,ind),reshape(u{j},[npv nc ng nj]));                            
        end
    end
end

