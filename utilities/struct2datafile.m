function struct2datafile(mesh,master,app,UDG,UH,filename)

% open binary file to write
fid = fopen(filename,'wb');

% mesh structure
fwritearray(fid, int64(mesh.nd));
fwritearray(fid, int64(mesh.porder));
fwritearray(fid, int64(mesh.elemtype));
fwritearray(fid, int64(mesh.nodetype));
fwritearray(fid, mesh.dgnodes);
fwritearray(fid, int64(mesh.bf));
fwritearray(fid, int64(mesh.elcon));
fwritearray(fid, int64(mesh.f));
fwritearray(fid, int64(mesh.t2f));

% master structure
fwritearray(fid, master.shapmv);
fwritearray(fid, master.shapvt);
fwritearray(fid, master.shapvg);
fwritearray(fid, master.shapvgdotshapvl);
fwritearray(fid, master.shapmf);
fwritearray(fid, master.shapft);
fwritearray(fid, master.shapfg);
fwritearray(fid, master.shapfgdotshapfc);
fwritearray(fid, int64(master.perm));
fwritearray(fid, int64(master.permgeom));

% solution structure
fwritearray(fid, UDG);
fwritearray(fid, UH);

% app structure
checkfield(app,'nc','app');
checkfield(app,'ncu','app');
checkfield(app,'ncq','app');
checkfield(app,'nch','app');
checkfield(app,'arg','app');
checkfield(app,'bcm','app');
checkfield(app,'bcs','app');
checkfield(app,'bcd','app');
checkfield(app,'bcv','app');
checkfield(app,'time','app');
checkfield(app,'dt','app');
checkfield(app,'fc_q','app');
checkfield(app,'fc_u','app');
checkfield(app,'fc_p','app');
checkfield(app,'tdep','app');
checkfield(app,'wave','app');
checkfield(app,'adjoint','app');
checkfield(app,'getdqdg','app');
checkfield(app,'appname','app');

fwritearray(fid, int64(app.nc));
fwritearray(fid, int64(app.ncu));
fwritearray(fid, int64(app.ncq));
fwritearray(fid, int64(app.nch));
fwritearray(fid, cell2mat(app.arg));
fwritearray(fid, int64(app.bcm));
fwritearray(fid, app.bcs);
fwritearray(fid, int64(app.bcd));
fwritearray(fid, app.bcv);
fwritearray(fid, app.time);
fwritearray(fid, app.dt);
fwritearray(fid, app.fc_q);
fwritearray(fid, app.fc_u);
fwritearray(fid, app.fc_p);
fwritearray(fid, logical(app.tdep));
fwritearray(fid, logical(app.wave));
fwritearray(fid, logical(app.adjoint));
fwritearray(fid, logical(app.getdqdg));
fwritearray(fid, app.appname);

% close the file
fclose(fid);

% check 
[mesh1,master1,app1,UDG1,UH1] = datafile2struct(filename);

fields = fieldnames(mesh1);
for i = 1:length(fields)         
    a1 = getfield(mesh1,fields{i});
    a = getfield(mesh,fields{i});    
    if max(abs(a(:)-a1(:))) > 1e-15
        error('something wrong');
    end
end

fields = fieldnames(master1);
for i = 1:length(fields)         
    a1 = getfield(master1,fields{i});
    a = getfield(master,fields{i});    
    if max(abs(a(:)-a1(:))) > 1e-15
        error('something wrong');
    end
end

if max(abs(UDG(:)-UDG1(:))) > 1e-15    
    error('something wrong');
end
if max(abs(UH(:)-UH1(:))) > 1e-15    
    error('something wrong');
end

fields = fieldnames(app1);
for i = 1:length(fields)         
    a1 = getfield(app1,fields{i});
    a = getfield(app,fields{i});   
    if iscell(a)
        a = cell2mat(a);
    end
    if max(abs(a(:)-a1(:))) > 1e-15
        error('something wrong');
    end    
end


function checkfield(mystruct,fieldstr,structstr)

if isfield(mystruct,fieldstr) == 0 || isempty(getfield(mystruct,fieldstr))
    str = ['Require nonempty field ' fieldstr ' in ' structstr ' structure'];
    error(str);
end

function fwritearray(fid, a)

sz  = int64(size(a));
dim = int64(length(sz));

if isa(a,'double')
    cln = 0;
elseif isa(a,'int64')
    cln = 1;
elseif isa(a,'logical')
    cln = 2;    
elseif isa(a,'char')
    cln = 3;        
else
    str = ['does not support this datatype' cl];
    error(str);
end
cln = int64(cln);

fwrite(fid,dim,'int64');
fwrite(fid,sz,'int64');
fwrite(fid,cln,'int64');
fwrite(fid,a,class(a));
