function dmd = readDMDSTRUCTfromBinaryFile(filename)

fileID = fopen([filename,'.bin'],'r');
N = fread(fileID,1,'int');
ndims = fread(fileID,N-1,'int');

% dmd structure
dmd.my_rank = ndims(1);
dmd.nbsd = fread(fileID,ndims(2),'int');
dmd.intelem = fread(fileID,ndims(3),'int');
dmd.intent = fread(fileID,ndims(4),'int');
dmd.elempart = fread(fileID,ndims(5),'int');
dmd.elempartpts = fread(fileID,ndims(6),'int');
dmd.entpart = fread(fileID,ndims(7),'int');
dmd.entpartpts = fread(fileID,ndims(8),'int');
dmd.elemrecv = fread(fileID,ndims(9),'int');
dmd.elemrecvpts = fread(fileID,ndims(10),'int');
dmd.elemsend = fread(fileID,ndims(11),'int');
dmd.elemsendpts = fread(fileID,ndims(12),'int');        
dmd.entrecv = fread(fileID,ndims(13),'int');
dmd.entrecvpts = fread(fileID,ndims(14),'int');
dmd.entsend = fread(fileID,ndims(15),'int');
dmd.entsendpts = fread(fileID,ndims(16),'int');
dmd.vecrecv = fread(fileID,ndims(17),'int');
dmd.vecrecvpts = fread(fileID,ndims(18),'int');
dmd.vecsend = fread(fileID,ndims(19),'int');
dmd.vecsendpts = fread(fileID,ndims(20),'int');
dmd.matrecv = fread(fileID,ndims(21),'int');
dmd.matrecvpts = fread(fileID,ndims(22),'int');
dmd.matsend = fread(fileID,ndims(23),'int');
dmd.matsendpts = fread(fileID,ndims(24),'int');
dmd.rowent2elem = fread(fileID,ndims(25),'int');
dmd.colent2elem = fread(fileID,ndims(26),'int');
dmd.rowent2ent = fread(fileID,ndims(27),'int');
dmd.colent2ent = fread(fileID,ndims(28),'int');
dmd.bcrs_rowent2elem = fread(fileID,ndims(29),'int');
dmd.bcrs_colent2elem = fread(fileID,ndims(30),'int');
dmd.bcrs_rowent2ent = fread(fileID,ndims(31),'int');
dmd.bcrs_colent2ent = fread(fileID,ndims(32),'int');        
dmd.ent2ind = fread(fileID,ndims(33),'int');
dmd.elcon = fread(fileID,ndims(34),'int');
dmd.t2f = fread(fileID,ndims(35),'int');
dmd.elemmap = fread(fileID,ndims(36),'int');
dmd.entmap = fread(fileID,ndims(37),'int');

fclose(fileID);
