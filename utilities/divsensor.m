function sdg = divsensor(mesh,mastersubgrid,udg,uhat)

if nargin>3
    qdg = getq(mesh,mastersubgrid,udg,uhat);  
end

for j = 1:length(udg)               
    if isempty(udg{j})==0             
        if nargin>3
            sdg{j} = sensor(cat(2,udg{j},qdg{j}));
        else
            sdg{j} = sensor(udg{j});
        end
    end        
end
sdg = udgproj(mesh,mastersubgrid,sdg);
