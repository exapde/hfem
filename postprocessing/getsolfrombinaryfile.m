function [UDG,UH] = getsolfrombinaryfile(filename,nproc,npv,nc,npf,nch,hybrid)

if nproc==1
    fileID = fopen([filename,'.bin'],'r');
    ndims = fread(fileID,2,'double');
    ne = ndims(1)/(npv*nc);
    UDG = reshape(fread(fileID,ndims(1),'double'),[npv nc ne]);    
    UH = reshape(fread(fileID,ndims(2),'double'),nch,[]);
    fclose(fileID);
else
%     fileID = fopen([filename,'elempart.bin'],'r');
%     elempart = fread(fileID,'int');
%     fclose(fileID);
%     nelem = elempart(1:nproc); elempart = elempart((nproc+1):end)+1;
%     n1 = cumsum([0; nelem]); 
%     
%     fileID = fopen([filename,'entpart.bin'],'r');
%     entpart = fread(fileID,'int');
%     fclose(fileID);
%     nent = entpart(1:nproc); entpart = entpart((nproc+1):end)+1;
%     n2 = cumsum([0; nent]); 
    
    elempart = []; entpart = [];
    nelem = zeros(nproc,1);
    nent = zeros(nproc,1);
    for i = 1:nproc
        fileID = fopen([filename 'elempart_np' num2str(i-1) '.bin'],'r');
        elemparti = fread(fileID,'int');
        nelem(i) = length(elemparti);
        elempart = [elempart; elemparti];
        fclose(fileID);
        
        fileID = fopen([filename 'entpart_np' num2str(i-1) '.bin'],'r');
        entparti = fread(fileID,'int');
        nent(i) = length(entparti);
        entpart = [entpart; entparti];
        fclose(fileID);
    end
    elempart = elempart + 1;
    entpart = entpart+1;
    
    n1 = cumsum([0; nelem]); 
    n2 = cumsum([0; nent]);     
    
    UDG = zeros(npv,nc,n1(end));
    if strcmp(hybrid,'hdg')==1           
        UH = zeros(nch,npf,n2(end));
    else
        UH = zeros(nch,n2(end));
    end        
    
    for i = 1:nproc        
        fileID = fopen([filename '_np' num2str(i-1) '.bin'],'r');
        ndims = fread(fileID,2,'double');                               
        
        ne = nelem(i);                        
        UDG(:,:,elempart((n1(i)+1):(n1(i)+ne))) = reshape(fread(fileID,ndims(1),'double'),[npv nc ne]);
        
        nf = nent(i);                
        if strcmp(hybrid,'hdg')==1                                    
            UH(:,:,entpart((n2(i)+1):(n2(i)+nf))) = reshape(fread(fileID,ndims(2),'double'),[nch npf nf]);            
        else            
            UH(:,entpart((n2(i)+1):(n2(i)+nf))) = reshape(fread(fileID,ndims(2),'double'),[nch nf]);            
        end
        
        fclose(fileID);
    end
end

