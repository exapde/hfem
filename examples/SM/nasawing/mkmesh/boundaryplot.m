function boundaryplot(p,t,pars)

nd = size(p,2);
npe = size(t,2);
elemtype = 0;
if npe==2^nd
    elemtype = 1;
end
    
f = mkt2f(t,elemtype);

patch('faces',f(:,1:end-2),'vertices',p,pars{:});                         

