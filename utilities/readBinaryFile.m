
function [master,mesh,app,UDG,UH] = readBinaryFile(filename)
% d: dimension
% e: element
% f: face
% v: vertex
% c: component
% p: polynomial
% g: gauss
% m: geometry
% t: time
% n: number 
% s: subdomain

fileID = fopen([filename,'.bin'],'r');

% dimensions
ndims = fread(fileID,100,'double');
mesh.nd = ndims(1);
mesh.ncd = ndims(2);
mesh.nfe = ndims(3);
mesh.nve = ndims(4);
mesh.nvf = ndims(5);
mesh.ne = ndims(6);
mesh.nf = ndims(7);
mesh.nv = ndims(8);
mesh.ndh = ndims(9); % degrees of UH
master.npe = ndims(10); % solution
master.npf = ndims(11);
master.nme = ndims(12); % geometry
master.nmf = ndims(13);
master.nge = ndims(14); % Gauss
master.ngf = ndims(15);
app.porder = ndims(16); % solution order
app.morder = ndims(17); % geometry order 
app.torder = ndims(18); % time order    
app.nstage = ndims(19); % number of DIRK stages
app.nc = ndims(20);
app.ncu = ndims(21);
app.ncq = ndims(22);
app.ncp = ndims(23);
app.nch = ndims(24);
app.ns  = ndims(25);
app.nb  = ndims(26);
app.ndt = ndims(27);
app.nparam = ndims(28);
app.nflag = ndims(29);
app.nfactor = ndims(30);

ne = mesh.ne;
nf = mesh.nf;
nd = mesh.nd;
ncd = mesh.ncd;
nfe = mesh.nfe;
ndh = mesh.ndh;
nvf = mesh.nvf;
npe = master.npe;
nge = master.nge;
nme = master.nme;
npf = master.npf;
ngf = master.ngf;
nmf = master.nmf;
nc = app.nc;
nch = app.nch;
%ncu = app.ncu;
ndt = app.ndt;
nb = app.nb;
nparam = app.nparam;
nflag = app.nflag;
nfactor = app.nfactor;

% solution
UDG = fread(fileID,npe*nc*ne,'double');
UH  = fread(fileID,nch*ndh,'double');

% master
master.plocvl = fread(fileID,npe*nd,'double');
master.gpvl = fread(fileID,nge*nd,'double');
master.gwvl = fread(fileID,nge,'double');
master.plocfc = fread(fileID,npf*(nd-1),'double');
master.gpfc = fread(fileID,ngf*(nd-1),'double');
master.gwfc = fread(fileID,ngf,'double');
master.shapvt = fread(fileID,nge*npe*(nd+1),'double');
master.shapvg = fread(fileID,npe*nge*(nd+1),'double');
master.shapvgdotshapvl = fread(fileID,npe*npe*nge*(nd+1),'double');
master.shapft = fread(fileID,ngf*npf*nd,'double');
master.shapfg = fread(fileID,npf*ngf*nd,'double');
master.shapfgdotshapfc = fread(fileID,npf*npf*ngf*nd,'double');
master.shapmv = fread(fileID,nge*nme*(nd+1),'double');
master.shapmf = fread(fileID,ngf*nmf*nd,'double');

% app
app.bcm = fread(fileID,nb,'double');
app.bcs = fread(fileID,nb*nch,'double');
app.bcd = fread(fileID,nb,'double');
app.bcv = fread(fileID,nb*nch,'double');
app.dt = fread(fileID,ndt,'double');
app.param = fread(fileID,nparam,'double');
app.flag = fread(fileID,nflag,'double');
app.factor = fread(fileID,nfactor,'double');

% mesh
mesh.dgnodes = fread(fileID,nme*ncd*ne,'double'); % local dgnodes on each subdomain
mesh.elcon = fread(fileID,npf*nfe*ne,'double');
mesh.bf = fread(fileID,nfe*ne,'double');
mesh.t2f = fread(fileID,ne*nfe,'double');
mesh.perm = fread(fileID,npf*nfe,'double');
mesh.permgeom = fread(fileID,nmf*nfe,'double');

fclose(fileID);

return;


