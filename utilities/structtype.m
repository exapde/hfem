function [mesh,master,app] = structtype(mesh,master,app)

% mesh structure
mesh.nd =  int64(mesh.nd);
mesh.porder = int64(mesh.porder);
mesh.elemtype = int64(mesh.elemtype);
mesh.nodetype = int64(mesh.nodetype);
mesh.bf = int64(mesh.bf);
mesh.elcon = int64(mesh.elcon);
mesh.f = int64(mesh.f);
mesh.t2f = int64(mesh.t2f);

% master structure
if min(master.perm(:))==1
    master.perm = master.perm - 1;
end
if min(master.permgeom(:))==1
    master.permgeom = master.permgeom-1;
end
master.perm = int64(master.perm);
master.permgeom = int64(master.permgeom);

app.nc = int64(app.nc);
app.ncu = int64(app.ncu);
app.ncq = int64(app.ncq);
app.nch = int64(app.nch);
app.bcm = int64(app.bcm);
app.bcd = int64(app.bcd);
app.tdep = logical(app.tdep);
app.wave = logical(app.wave);
app.adjoint = logical(app.adjoint);
app.getdqdg = logical(app.getdqdg);
%app.appname = char(app.appname);

