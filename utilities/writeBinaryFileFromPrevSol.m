
function writeBinaryFileFromPrevSol(structsFile,solFile,inputsFile,meshpFile)

% structsFile: File containing mesh, app, master, UDG, UH structs
% solFile: Binary file with previous solution. Should include time step
% inputsFile: Binary files to be read by C++ code
% meshpFile (optional. Only for parallel computations): File containing meshp structure

% Load structures:
load(structsFile);
if nargin > 3
    load(meshpFile);
    nproc = length(meshp);
else
    warning('meshpFile file was not provided. Serial computation is assumed.');
    nproc = 1;
end

if(app.wave==1); error('This function does not work for wave problems.'); end
if length(app.dt) ~= app.ndt; error('app.ndt value does not match app.dt length.'); end

% Read solution file
if strcmp(mesh.hybrid,'hdg'); UH = reshape(UH, [app.nch size(mesh.perm,1) mesh.nf]); end
if nproc>1
    for i = 1:nproc
        filename = [solFile, '_np' num2str(i-1)];
        fileID = fopen([filename,'.bin'],'r');
        data = fread(fileID,'double');    
        fclose(fileID);

        nei = sum(meshp{i}.elempartpts(1:2));
        nfi = sum(meshp{i}.entpartpts(1:2)); 
        inde = meshp{i}.elempart(1:nei);               
        indf = meshp{i}.entpart(1:nfi);

        npf = size(mesh.perm,1);
        nc = app.nc;    
        nch = app.nch;    
        npe = size(mesh.plocal,1);        
        n1 = npe*nc*nei;        
        
        UDG(:,:,inde) = reshape(data(1:n1),[npe nc nei]);
        if strcmp(mesh.hybrid,'hdg');
            UH(:,:,indf) = reshape(data(n1+1:end),[nch npf nfi]);
        else
            UH(:,indf) = reshape(data(n1+1:end),[nch nfi]);
        end
    end
else
    fileID = fopen([solFile,'.bin'],'r');
    data = fread(fileID,'double');    
    fclose(fileID);

    npf = size(mesh.perm,1);
    nc = app.nc;    
    nch = app.nch;    
    npe = size(mesh.plocal,1);    
    nf = mesh.nf;
    ne = mesh.ne;
    numEnt = mesh.cbsr_nrows;
    n1 = npe*nc*ne;        

    UDG(:,:,:) = reshape(data(1:n1),[npe nc ne]);
    if strcmp(mesh.hybrid,'hdg');
        UH(:,:,:) = reshape(data(n1+1:end),[nch npf nf]);
    else
        UH(:,:) = reshape(data(n1+1:end),[nch numEnt]);
    end
end


% Write binary files:
if nproc == 1
    meshp = [];
    
    ndims = zeros(100,1);
    ndims(1) = mesh.nd;
    ndims(2) = mesh.ncd; % components of dgnodes
    ndims(3) = mesh.nfe; % faces per element
    ndims(4) = mesh.nve; % vertices per element
    ndims(5) = mesh.nvf; % vertices per face
    ndims(6) = mesh.ne;
    ndims(7) = mesh.nf;
    ndims(8) = mesh.nv;
    ndims(9) = mesh.ndh; % degrees of UH
    ndims(10) = master.npe; % solution
    ndims(11) = master.npf;
    ndims(12) = master.nme; % geometry
    ndims(13) = master.nmf;
    ndims(14) = master.nge; % Gauss
    ndims(15) = master.ngf;
    ndims(16) = app.porder; % solution order
    ndims(17) = app.morder; % geometry order 
    ndims(18) = app.torder; % time order    
    ndims(19) = app.nstage; % number of DIRK stages
    ndims(20) = app.nc;
    ndims(21) = app.ncu;
    ndims(22) = app.ncq;
    ndims(23) = app.ncp;
    ndims(24) = app.nch;
    ndims(25) = app.ns; % number of subdomains
    ndims(26) = app.nb; % number of boundaries
    ndims(27) = app.ndt;
    ndims(28) = app.nparam;
    ndims(29) = app.nflag;
    ndims(30) = app.nfactor;
    ndims(31) = app.nproblem;
    ndims(32) = mesh.cbsr_nrows;
    ndims(33) = mesh.cbsr_nblks;
    ndims(34) = mesh.maxBlocksPerRow;
    ndims(35) = mesh.blkSize;
    
    ndims(40) = nproc;

    fileID = fopen([inputsFile,'.bin'],'w');

    % dimensions
    fwrite(fileID,ndims,'double');
    
    % master
    fwrite(fileID,master.plocvl,'double');
    fwrite(fileID,master.gpvl,'double');
    fwrite(fileID,master.gwvl,'double');
    fwrite(fileID,master.plocfc,'double');
    fwrite(fileID,master.gpfc,'double');
    fwrite(fileID,master.gwfc,'double');
    fwrite(fileID,master.shapvt,'double');
    fwrite(fileID,master.shapvg,'double');
    fwrite(fileID,master.shapvgdotshapvl,'double');
    fwrite(fileID,master.shapft,'double');
    fwrite(fileID,master.shapfg,'double');
    fwrite(fileID,master.shapfgdotshapfc,'double');
    fwrite(fileID,master.shapmv,'double');
    fwrite(fileID,master.shapmf,'double');

    % app    
    fwrite(fileID,app.bcm(:),'double');
    fwrite(fileID,app.bcs','double');
    fwrite(fileID,app.bcd(:),'double');
    fwrite(fileID,app.bcv','double');
    fwrite(fileID,app.dt(:),'double');
    fwrite(fileID,app.param(:),'double');
    fwrite(fileID,app.flag(:),'double');
    fwrite(fileID,app.factor(:),'double');
    fwrite(fileID,app.problem(:),'double');
    
    % mesh
    fwrite(fileID,mesh.dgnodes,'double'); % local dgnodes on each subdomain
    fwrite(fileID,mesh.elcon,'double');
    fwrite(fileID,mesh.bf,'double');
    fwrite(fileID,mesh.t2f,'double');
    fwrite(fileID,mesh.perm,'double');
    fwrite(fileID,mesh.permgeom,'double');
    
    % solution
    fwrite(fileID,UDG,'double');
    fwrite(fileID,UH,'double');    
    if app.wave==1
        fwrite(fileID,PDG,'double');
    end
    
    % sys
    fwrite(fileID,mesh.cbsr_rowpts,'double');
    fwrite(fileID,mesh.cbsr_colind,'double');

    fclose(fileID);
    
    % Compute memory requirements:
    ncu = app.ncu;
    nch = app.nch;
    ne = mesh.ne;
    nfe = mesh.nfe;
    npv = master.npe;
    npf = master.npf;
    blkSize = mesh.blkSize;
    cbsr_nblks = mesh.cbsr_nblks;

    memory = zeros(4,1);
    memory(1) = 8*blkSize*blkSize*cbsr_nblks/(1024^2);   % Hg
    memory(2) = 8*blkSize*blkSize*cbsr_nblks/(1024^2);   % Mg
    memory(3) = 8*blkSize*cbsr_nblks/(1024^2);             % v
    memory(4) = 8*npv*ncu*npf*nfe*nch*ne/(1024^2);       % DinvF
    memory(5) = 4*npv*ncu*npv*ncu*ne/(1024^2);       % Dinv
    memory(6) = 4*nch*npf*nfe*npv*ncu*ne/(1024^2);       % K
    
    disp(' ');
    disp(' ');
    disp('************************************************************');
    disp('************** Summary of memory requirements **************');
    disp('************************************************************');
    disp(' ');
    disp(['Hg: ', num2str(memory(1)), ' MB']);
    disp(['Mg: ', num2str(memory(2)), ' MB']);
    disp(['v: Restart x ', num2str(memory(3)),' MB']);
    disp(['DinvF: ', num2str(memory(4)),' MB']);
    disp(['Dinv: ', num2str(memory(5)),' MB']);
    disp(['K: ', num2str(memory(6)),' MB']);
    disp(' ');
    totalMemory1 = (memory(1) + memory(2) + memory(4)) / 1024;
    totalMemory2 = memory(3) / 1024;
    totalMemory3 = (memory(1) + memory(2) + memory(4) + memory(5) + memory(6)) / 1024;
    disp(['TOTAL MEMORY REQUIREMENTS (Newton): ', num2str(totalMemory1), ' GB  + Restart x ', num2str(totalMemory2), ' GB.']);
    disp(['TOTAL MEMORY REQUIREMENTS (quasi-Newton): ', num2str(totalMemory3), ' GB  + Restart x ', num2str(totalMemory2), ' GB.']);
    
elseif nproc > 1
    nelem = zeros(1,nproc);
    memory = zeros(nproc,6);
    for i = 1:nproc
        nelem(i) = meshp{i}.ne;
        fileID = fopen([inputsFile num2str(i-1) '.bin'],'w');
        
        ndims = zeros(100,1);
        ndims(1) = meshp{i}.nd;
        ndims(2) = meshp{i}.ncd; % components of dgnodes
        ndims(3) = meshp{i}.nfe; % faces per element
        ndims(4) = meshp{i}.nve; % vertices per element
        ndims(5) = meshp{i}.nvf; % vertices per face
        ndims(6) = meshp{i}.ne;
        ndims(7) = meshp{i}.nf;
        ndims(8) = meshp{i}.nv;
        ndims(9) = meshp{i}.ndh; % degrees of UH
        ndims(10) = master.npe; % solution
        ndims(11) = master.npf;
        ndims(12) = master.nme; % geometry
        ndims(13) = master.nmf;
        ndims(14) = master.nge; % Gauss
        ndims(15) = master.ngf;
        ndims(16) = app.porder; % solution order
        ndims(17) = app.morder; % geometry order 
        ndims(18) = app.torder; % time order    
        ndims(19) = app.nstage; % number of DIRK stages
        ndims(20) = app.nc;
        ndims(21) = app.ncu;
        ndims(22) = app.ncq;
        ndims(23) = app.ncp;
        ndims(24) = app.nch;
        ndims(25) = app.ns; % number of subdomains
        ndims(26) = app.nb; % number of boundaries
        ndims(27) = app.ndt;
        ndims(28) = app.nparam;
        ndims(29) = app.nflag;
        ndims(30) = app.nfactor;
        ndims(31) = app.nproblem;
        ndims(32) = meshp{i}.cbsr_nrows;
        ndims(33) = meshp{i}.cbsr_nblks;
        ndims(34) = meshp{i}.maxBlocksPerRow;
        ndims(35) = meshp{i}.blkSize;
        ndims(36) = meshp{i}.BJ_nrows;
        ndims(37) = meshp{i}.BJ_nblks;
        ndims(38) = meshp{i}.BK_nrows;
        ndims(39) = meshp{i}.BK_nblks;
        ndims(40) = nproc;
        ndims(41) = length(meshp{i}.entpartpts(:));
        ndims(42) = length(meshp{i}.entrecv(:));
        ndims(43) = length(meshp{i}.entsend(:));
        ndims(44) = length(meshp{i}.elempartpts(:));
        ndims(45) = length(meshp{i}.elemrecv(:));
        ndims(46) = length(meshp{i}.elemsend(:));
        ndims(47) = length(meshp{i}.nbsd(:));
        ndims(48) = length(mesh.cbsr_rowpts(:));    % This variable is no longer required in the C++ cpde
        ndims(49) = length(mesh.cbsr_colind(:));    % This variable is no longer required in the C++ code
        ndims(50) = length(meshp{i}.matrecv(:));
        ndims(51) = length(meshp{i}.matsend(:));
        
        % dimensions
        fwrite(fileID,ndims,'double');
        N = 100;
        
        % master
        fwrite(fileID,master.plocvl,'double');
        fwrite(fileID,master.gpvl,'double');
        fwrite(fileID,master.gwvl,'double');
        fwrite(fileID,master.plocfc,'double');
        fwrite(fileID,master.gpfc,'double');
        fwrite(fileID,master.gwfc,'double');
        fwrite(fileID,master.shapvt,'double');
        fwrite(fileID,master.shapvg,'double');
        fwrite(fileID,master.shapvgdotshapvl,'double');
        fwrite(fileID,master.shapft,'double');
        fwrite(fileID,master.shapfg,'double');
        fwrite(fileID,master.shapfgdotshapfc,'double');
        fwrite(fileID,master.shapmv,'double');
        fwrite(fileID,master.shapmf,'double');
        
        N = N+length(master.plocvl(:))+length(master.gpvl(:))+length(master.gwvl(:))...
            +length(master.plocfc(:))+length(master.gpfc(:))+length(master.gwfc(:))...
            +length(master.shapvt(:))+length(master.shapvg(:))+length(master.shapvgdotshapvl(:))...
            +length(master.shapft(:))+length(master.shapfg(:))+length(master.shapfgdotshapfc(:))...
            +length(master.shapmv(:))+length(master.shapmf(:));
        
        % app        
        fwrite(fileID,app.bcm(:),'double');
        fwrite(fileID,app.bcs','double');
        fwrite(fileID,app.bcd(:),'double');
        fwrite(fileID,app.bcv','double');
        fwrite(fileID,app.dt(:),'double');
        fwrite(fileID,app.param(:),'double');
        fwrite(fileID,app.flag(:),'double');
        fwrite(fileID,app.factor(:),'double');
        fwrite(fileID,app.problem(:),'double');

        N = N+length(app.bcm(:))+length(app.bcs(:))+length(app.bcd(:))+length(app.bcv(:))+length(app.dt(:))+...
            length(app.param(:))+length(app.flag(:))+length(app.factor(:))+length(app.problem(:));
        
        % mesh
        fwrite(fileID,meshp{i}.dgnodes,'double'); % local dgnodes on each subdomain
        fwrite(fileID,meshp{i}.elcon,'double');
        fwrite(fileID,meshp{i}.bf,'double');
        fwrite(fileID,meshp{i}.t2f,'double');
        fwrite(fileID,meshp{i}.perm,'double');
        fwrite(fileID,meshp{i}.permgeom,'double');
        
        % solution
        fwrite(fileID,UDG(:,:,meshp{i}.elempart),'double');   
        N = N+numel(UDG(:,:,meshp{i}.elempart));
        if strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg') || strcmp(mesh.hybrid,'hdg')
            fwrite(fileID,UH(:,meshp{i}.entpart),'double');
            N = N + numel(UH(:,meshp{i}.entpart));
        elseif strcmp(mesh.hybrid,'hdg')
            UH = reshape(UH,[app.ncu master.npf mesh.nf]);
            fwrite(fileID,UH(:,:,meshp{i}.entpart),'double');
            N = N + numel(UH(:,:,meshp{i}.entpart));
        end
        if app.wave==1
            fwrite(fileID,PDG(:,:,meshp{i}.elempart),'double');
            N = N + numel(PDG(:,:,meshp{i}.elempart));
        end
        
        % sys
        fwrite(fileID,meshp{i}.cbsr_rowpts,'double');
        fwrite(fileID,meshp{i}.cbsr_colind,'double');
%         if i == 1
%              fwrite(fileID,mesh.globalEnt2entStart,'double');
%              fwrite(fileID,mesh.globalEnt2ent,'double');
%         end
        fwrite(fileID,meshp{i}.BJ_rowpts,'double');
        fwrite(fileID,meshp{i}.BJ_colind,'double');
        fwrite(fileID,meshp{i}.BK_rowpts,'double');
        fwrite(fileID,meshp{i}.BK_colind,'double');
        fwrite(fileID,meshp{i}.BK_rowind,'double');
        fwrite(fileID,meshp{i}.entpart,'double');
        fwrite(fileID,meshp{i}.entpartpts,'double');
        fwrite(fileID,meshp{i}.entrecv,'double');
        fwrite(fileID,meshp{i}.entrecvpts,'double');
        fwrite(fileID,meshp{i}.entsend,'double');
        fwrite(fileID,meshp{i}.entsendpts,'double');
        fwrite(fileID,meshp{i}.elempart,'double');
        fwrite(fileID,meshp{i}.elempartpts,'double');
        fwrite(fileID,meshp{i}.elemrecv,'double');
        fwrite(fileID,meshp{i}.elemrecvpts,'double');
        fwrite(fileID,meshp{i}.elemsend,'double');
        fwrite(fileID,meshp{i}.elemsendpts,'double');
        fwrite(fileID,meshp{i}.nbsd,'double'); 
        fwrite(fileID,meshp{i}.matrecv,'double');
        fwrite(fileID,meshp{i}.matrecvpts,'double');
        fwrite(fileID,meshp{i}.matsend,'double');
        fwrite(fileID,meshp{i}.matsendpts,'double');
        
%         if i == 1
%             for j=1:nproc
%              fwrite(fileID,meshp{j}.ent2entWeightLen,'double');
%             end
%         end
    
        N=N+length(meshp{i}.dgnodes(:))+length(meshp{i}.elcon(:))+length(meshp{i}.bf(:))+...
            length(meshp{i}.t2f(:))+length(meshp{i}.perm(:))+length(meshp{i}.permgeom(:))+...
            length(meshp{i}.cbsr_rowpts(:))+length(meshp{i}.cbsr_colind(:))+...
            length(meshp{i}.BJ_rowpts(:))+length(meshp{i}.BJ_colind(:))+length(meshp{i}.BK_rowpts(:))+...
            length(meshp{i}.BK_colind(:))+length(meshp{i}.BK_rowind(:))+...
            length(meshp{i}.entpart(:))+length(meshp{i}.entpartpts(:))+...
            length(meshp{i}.elempart(:))+length(meshp{i}.elempartpts(:))+...
            length(meshp{i}.entrecv(:))+length(meshp{i}.entrecvpts(:))+...
            length(meshp{i}.elemrecv(:))+length(meshp{i}.elemrecvpts(:))+...
            length(meshp{i}.entsend(:))+length(meshp{i}.entsendpts(:))+...
            length(meshp{i}.elemsend(:))+length(meshp{i}.elemsendpts(:))+length(meshp{i}.nbsd(:))+...
            length(meshp{i}.matrecv(:))+length(meshp{i}.matrecvpts(:))+...
            length(meshp{i}.matsend(:))+length(meshp{i}.matsendpts(:));
        
%         if i == 1
%             N = N + length(mesh.globalEnt2entStart) + length(mesh.globalEnt2ent) + nproc;
%         end

        % Compute memory requirements:
        ncu = app.ncu;
        nch = app.nch;
        ne = meshp{i}.ne;
        nfe = meshp{i}.nfe;
        npv = master.npe;
        npf = master.npf;
        blkSize = meshp{i}.blkSize;
        cbsr_nblks = meshp{i}.cbsr_nblks;
        
        memory(i,1) = 8*blkSize*blkSize*cbsr_nblks/(1024^2);   % Hg
        memory(i,2) = 8*blkSize*blkSize*cbsr_nblks/(1024^2);   % Mg
        memory(i,3) = 8*blkSize*cbsr_nblks/(1024^2);             % v
        memory(i,4) = 8*npv*ncu*npf*nfe*nch*ne/(1024^2);       % DinvF
        memory(i,5) = 4*npv*ncu*npv*ncu*ne/(1024^2);       % Dinv
        memory(i,6) = 4*nch*npf*nfe*npv*ncu*ne/(1024^2);       % K
        
        fclose(fileID);
    end
    
    disp(' ');
    disp(' ');
    disp('************************************************************');
    disp('************** Summary of memory requirements **************');
    disp('************************************************************');
    disp(' ');
    for i=1:nproc
        disp(['Proc. #',num2str(i),' Hg: ',num2str(memory(i,1)), ...
            ' MB   ||   Mg: ',num2str(memory(i,2)),' MB   ||   v: Restart x ', ...
            num2str(memory(i,3)),' MB   ||   DinvF: ',num2str(memory(i,4)),' MB   ||   Dinv: ',...
            num2str(memory(i,5)),' MB   ||   K: ',num2str(memory(i,6)),' MB']);
    end
    disp(' ');
    totalMemory1 = (sum(memory(:,1)) + sum(memory(:,2)) + sum(memory(:,4))) / 1024;
    totalMemory2 = sum(memory(:,3)) / 1024;
    totalMemory3 = (sum(memory(:,1)) + sum(memory(:,2)) + sum(memory(:,4)) + sum(memory(:,5)) + sum(memory(:,6))) / 1024;
    disp(['TOTAL MEMORY REQUIREMENTS (Newton): ', num2str(totalMemory1), ' GB + Restart x ', num2str(totalMemory2), ' GB.']);
    disp(['TOTAL MEMORY REQUIREMENTS (quasi-Newton): ', num2str(totalMemory3), ' GB + Restart x ', num2str(totalMemory2), ' GB.']);
else
    error('Invalid number of processors');
end

% Save structures into .mat file:
save([inputsFile, '.mat']);
