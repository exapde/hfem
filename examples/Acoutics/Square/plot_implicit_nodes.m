figure()
meshplot(mesh, [0 0 0])
x = reshape(mesh.dgnodes(:,1,ind), [impE*10,1]);
y = reshape(mesh.dgnodes(:,2,ind), [impE*10,1]);
hold on
plot(x,y,'*')

% % IMEX Boundary
% x2 = reshape(mesh.dgnodes(:,1,indf), [faceimex*10,1]);
% y2 = reshape(mesh.dgnodes(:,2,indf), [faceimex*10,1]);
% plot(x2,y2,'*')