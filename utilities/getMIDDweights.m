
function globalEnt2entWeight = getMIDDweights(filename,mesh, meshp)

nproc = length(meshp);
globalEnt2entWeight = zeros(length(mesh.cbsr_colind),1);
for i=1:nproc
    fileID = fopen([filename,'_ent2entWeight_np',num2str(i-1),'.bin'],'r');
    meshp{i}.ent2entWeight = fread(fileID,meshp{i}.BJ_nblks+meshp{i}.BK_nblks,'double');
    fclose(fileID);
    
    rowIndices = zeros(meshp{i}.BJ_nblks+meshp{i}.BK_nblks,1);
    startt = 1;
    for j=1:(meshp{i}.entpartpts(1)+meshp{i}.entpartpts(2))
        i1 = mesh.cbsr_rowpts(meshp{i}.entpart(j))+1;
        i2 = mesh.cbsr_rowpts(meshp{i}.entpart(j)+1);
        rowIndices(startt:(startt+i2-i1)) = i1:i2;
        startt = startt + i2 - i1 + 1;
    end
    if startt ~= (meshp{i}.BJ_nblks+meshp{i}.BK_nblks+1); error('Something wrong.'); end
    if min(rowIndices) < 1; error('Something wrong.'); end

    globalEnt2entWeight(rowIndices) = meshp{i}.ent2entWeight(:);
end

end