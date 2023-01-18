% clear all
load example2dMesh.mat
nlayers = 4;
thickness = 1;

[mesh3d.p, mesh3d.t] = QtoT(mesh2d.p, mesh2d.t, thickness, nlayers);

mesh3d.f = mkt2f(mesh3d.t);
figure(2)
meshplot(mesh3d)

