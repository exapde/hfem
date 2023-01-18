
function reduceSizeBinaryFile(mesh, app, dmd, inName, outName, nproc, tInit, tEnd, tFreq, removeOldFiles)

% Takes a binary file with (u,q,uh) and creates another with (u,uh) only

% inName = 'out_epp387Re100kAoA4d2Medium_40c';
% outName = '2out_epp387Re100kAoA4d2Medium_40c';
% nproc = 128;
% tInit = 1100;
% tEnd = 1100;
% tFreq = 10;

if nargin < 10; removeOldFiles = 0; end

ncu = app.ncu;
nc = app.nc;
npe = size(mesh.plocal,1);

% Loop over all binary files:
for n=tInit:tFreq:tEnd
    disp(['Time-step No. ',num2str(n)]);
    for i = 1:nproc
        % Read binary file with u and q:
        filename = [inName,'_t' num2str(n) '_np' num2str(i-1)];
        fileID = fopen([filename,'.bin'],'r');
        
        dataWithQ = fread(fileID,'double');    
        fclose(fileID);
        
        nei = sum(dmd{i}.elempartpts(1:2));             
        n1 = npe*nc*nei;        
        data_tmp = reshape(dataWithQ(1:n1),[npe nc nei]);
        
        % Write binary file with u only
        filename = [outName,'_t' num2str(n) '_np' num2str(i-1)];
        fileID = fopen([filename,'.bin'],'w');
        
        dataNoQ = [reshape(data_tmp(:,1:ncu,:),[],1); dataWithQ(n1+1:end)];
            
        fwrite(fileID,dataNoQ,'double');
        fclose(fileID);
        
        % Delete old file (if required)
        if removeOldFiles == 1;
            filename = [inName,'_t' num2str(n) '_np' num2str(i-1) '.bin'];
            delete(filename);
        end
    end
end
