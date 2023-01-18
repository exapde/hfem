function [pb,tb,eb,pp] = mkbeam2(p, pp, tp, ts, tn, ind, tol, isplot)

[beamline,normal,tangent] = getbeamline(p, ts, ind, tol);

[np, nd, ns] = size(normal);

n = size(beamline,1);
joinp = zeros(nd, 4, np);
joint = zeros(np, 4);
edge = zeros(nd, 2, ns, np);
for i = 1:n            
    [pi, joint(i,:), ei] = mkjoint(beamline(i,1:2)', reshape(tangent(i,1:2,:),2,[]), tn(ind));    
    joinp(1:2,:,i) = pi; 
    joinp(3,:,i) = beamline(i,3);
    joint(i,:) = joint(i,:) + (i-1)*4;    
    edge(1:2,:,:,i) = ei;
    edge(3,:,:,i) = beamline(i,3);
end

pb = reshape(permute(joinp, [2 3 1]),[4*np nd]);
tb = zeros(np-1,8);
for i = 1:np-1
    tb(i,:) = [joint(i,:) joint(i+1,:)];
end
eb = reshape(permute(edge, [2 4 1 3]),[2*np nd ns]);

for i = 1:ns
    pp{ind(i)} = fixp(pp{ind(i)}, eb(:,:,i));
end   

if nargin<8
    isplot=0;
end

if isplot==1
figure(2); clf;
hold on;
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},{'facecolor','g','edgecolor','k','Linew',1});     
end
for i = 1:length(ind)   
    plot3(eb(:,1,i),eb(:,2,i),eb(:,3,i),'ob');
end
boundaryplot(pb,tb,{'facecolor','r','edgecolor','k','Linew',1}); 
end



