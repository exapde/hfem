function [mesh,master,app,UDG,UH] = datafile2struct(filename)

% open binary file to write
fid = fopen(filename,'r');

% mesh structure
mesh.nd = freadarray(fid);
mesh.porder = freadarray(fid);
mesh.elemtype = freadarray(fid);
mesh.nodetype = freadarray(fid);
mesh.dgnodes = freadarray(fid);
mesh.bf = freadarray(fid);
mesh.elcon = freadarray(fid);
mesh.f = freadarray(fid);
mesh.t2f = freadarray(fid);

% master structure
master.shapmv = freadarray(fid);
master.shapvt = freadarray(fid);
master.shapvg = freadarray(fid);
master.shapvgdotshapvl = freadarray(fid);
master.shapmf = freadarray(fid);
master.shapft = freadarray(fid);
master.shapfg = freadarray(fid);
master.shapfgdotshapfc = freadarray(fid);
master.perm = freadarray(fid);
master.permgeom = freadarray(fid);

% solution structure
UDG = freadarray(fid);
UH = freadarray(fid);

% app structure
app.nc = freadarray(fid);
app.ncu = freadarray(fid);
app.ncq = freadarray(fid);
app.nch = freadarray(fid);
app.arg = freadarray(fid);
app.bcm = freadarray(fid);
app.bcs = freadarray(fid);
app.bcd = freadarray(fid);
app.bcv = freadarray(fid);
app.time = freadarray(fid);
app.dt = freadarray(fid);
app.fc_q = freadarray(fid);
app.fc_u = freadarray(fid);
app.fc_p = freadarray(fid);
app.tdep = freadarray(fid);
app.wave = freadarray(fid);
app.adjoint = freadarray(fid);
app.getdqdg = freadarray(fid);
app.appname = freadarray(fid);

% close the file
fclose(fid);


function [a,dim,sz,cln] = freadarray(fid)

dim = fread(fid,1,'int64');
sz = fread(fid,dim,'int64');
cln = fread(fid,1,'int64');
sz = reshape(sz,1,dim);

if cln == 0 % double
    cl = class(1.0);
elseif cln == 1 % 64-bit integer
    cl = class(int64(1));
elseif cln == 2 % boolean
    cl = class(false);
elseif cln == 3 % char
    cl = class('a');    
else
    str = 'does not support this datatype';
    error(str);
end

a = fread(fid,prod(sz),cl);
a = reshape(a,sz);

