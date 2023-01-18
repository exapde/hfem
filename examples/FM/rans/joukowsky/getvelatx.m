function v = getvelatx(udg,pdg,pm,x)

udg = permute(udg,[1 3 2]);
pdg = permute(pdg,[1 3 2]);

nx = size(x,1);
v = 0*x;
for i = 1:nx
    d = (pm(:,1)-x(i,1)).^2 + (pm(:,2)-x(i,2)).^2;
    dmin = min(d);
    ind = find(d<20*dmin);     
    
    in = 0;
    for j = 1:length(ind)
        k = ind(j);        
        in = PointInTriangle(x(i,:), pdg(1,:,k), pdg(2,:,k), pdg(3,:,k));
        if in==1            
            v(i,:) = (udg(1,:,k)+udg(2,:,k)+udg(3,:,k))/3;    
%             figure(1); clf;
%             plot(pdg(:,1,k),pdg(:,2,k),'ob');
%             hold on; 
%             plot(x(i,1),x(i,2),'xr');
            break;
        end
    end
    
    if in==0
        error('cannot find');
    end
end





