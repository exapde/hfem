function [newsg,a,q,u,s] = newssubgrid(mesh,mastersubgrid,udg,uh,c1,nlayer)

ns = length(mastersubgrid);

u = udgproj(mesh,mastersubgrid,udg);
if exist('uh','var')
    ncu = size(uh,1);    
    q = getq(mesh,mastersubgrid{1},u,reshape(uh(:,mesh.elcon),ncu,[],size(u,3)));
    s = sensor(cat(2,u,q));    
else
    s = sensor(u);    
end
a = mean(squeeze(s));

ne    = length(mesh.subgrids);
newsg = ones(ne,1);

% update new subgrids
ind0 = find(a <= c1); % p0 elements
if ns==2
    newsg(ind0) = 2; 
else
    ind2 = find(-0.1<=a & a<=10);               % high-order elements
    ind1 = setdiff((1:ne)',union(ind0,ind2)); % p1 elements    
    newsg(ind0) = 3;
    newsg(ind1) = 2;
end

t2t = mkt2t(mesh.t,mesh.elemtype);
for jj=1:nlayer    
    ind = find(newsg==ns);
    for i=1:length(ind)    
        t = t2t(ind(i),:);
        k = t > 0;    
        m = t(k);
        j = ismember(m,ind);
        newsg(m(j==0))=ns;    
    end
end

