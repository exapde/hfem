
nfe = size(master.perm,2);
nch = 4;
ncu = 4;
nc = 12;
ncd = size(mesh.dgnodes,2);
paramLength = 3;
param = {1,2,3};
flagsLength = 1;
flags = [1];
factorsLength = 2;
factors = [1,2];
numBoundaries = 2;
UDG = zeros(master.npv,nc,mesh.ne);
UH = zeros(master.npf,ncu,mesh.nf);
SH = zeros(master.npv,nc,mesh.ne);
bcm = [2,1];
ui = [1,1,0,20];
bcs = [ui,ui];
time = 1.0;

%%%%%%
%%%%%%

fileName = 'test';

fileID = fopen([fileName,'.bin'],'w');

data = [mesh.nd;mesh.ne;mesh.nf;master.npv;master.npf;nfe;master.ngv;master.ngf;ncu;nc;nch;ncd;paramLength;flagsLength;factorsLength;numBoundaries];
data = [data; reshape(UDG,[],1)];   % npv x nc x ne
data = [data; reshape(UH,[],1)];    % npf x ncu x nf
data = [data; reshape(SH,[],1)];    % npv x nc x ne
data = [data; reshape(mesh.dgnodes,[],1)];  % npv x ncd x ne
data = [data; reshape(mesh.bf,[],1)];   % nfe? x ne
data = [data; reshape(master.shapvt,[],1)]; % ngv x npv x nd+1
data = [data; reshape(master.shapvg,[],1)]; % npv x ngv x nd+1
data = [data; reshape(master.shapvgdotshapvl,[],1)];    % npv*npv x ngv x nd+1
data = [data; reshape(master.shapft,[],1)]; % ngf x npf x nd
data = [data; reshape(master.shapfg,[],1)]; % npf x ngf x nd
data = [data; reshape(master.shapfgdotshapfc,[],1)];    % npf*npf x ngf x nd
data = [data; reshape(master.perm,[],1)];       % npf x nfe
data = [data; reshape(cell2mat(param),[],1)];    % paramLength
data = [data; reshape(flags,[],1)];                 % flagsLength
data = [data; reshape(factors,[],1)];           % factorslength
data = [data; reshape(bcm,[],1)];           % numBoundaries
data = [data; reshape(bcs,[],1)];           % ncu x numBoundaries
data = [data; time];                        % 1
    
fileID = fopen([fileName,'.bin'],'w');
fwrite(fileID,data,'double');
fclose(fileID);

%%%%%
%%%%%

fileID = fopen([fileName,'.bin'],'r');
data = fread(fileID,'double');
fclose(fileID);

nd = data(1);
ne = data(2);
nf = data(3);
npv = data(4);
npf = data(5);
nfe = data(6);
ngv = data(7);
ngf = data(8);
ncu = data(9);
nc = data(10);
nch = data(11);
ncd = data(12);
paramLength = data(13);
flagsLength = data(14);
factorsLength = data(15);
numBoundaries = data(16);

iInit = 17;
iEnd = iInit + npv*nc*ne - 1;
UDGread = reshape(data(iInit:iEnd),[npv,nc,ne]);
iInit = iEnd + 1;
iEnd = iInit + npf*nch*nf - 1;
UHread = reshape(data(iInit:iEnd),[npf,nch,nf]);
iInit = iEnd + 1;
iEnd = iInit + npv*nc*ne - 1;
SHread = reshape(data(iInit:iEnd),[npv,nc,ne]);

iInit = iEnd + 1;
iEnd = iInit + npv*ncd*ne - 1;
dgnodesRead = reshape(data(iInit:iEnd),[npv,ncd,ne]);

iInit = iEnd + 1;
iEnd = iInit + nfe*ne - 1;
bfRead = reshape(data(iInit:iEnd),[nfe,ne]);

iInit = iEnd + 1;
iEnd = iInit + (ngv*npv*(nd+1)) - 1;
shapvtRead = reshape(data(iInit:iEnd),[ngv,npv,nd+1]);

iInit = iEnd + 1;
iEnd = iInit + npv*ngv*(nd+1) - 1;
shapvgRead = reshape(data(iInit:iEnd),[npv,ngv,nd+1]);

iInit = iEnd + 1;
iEnd = iInit + npv*npv*ngv*(nd+1) - 1;
shapvgdotshapvlRead = reshape(data(iInit:iEnd),[npv*npv,ngv,nd+1]);

iInit = iEnd + 1;
iEnd = iInit + ngf*npf*nd - 1;
shapftRead = reshape(data(iInit:iEnd),[ngf,npf,nd]);

iInit = iEnd + 1;
iEnd = iInit + npf*ngf*nd - 1;
shapfgRead = reshape(data(iInit:iEnd),[npf,ngf,nd]);

iInit = iEnd + 1;
iEnd = iInit + npf*npf*ngf*nd - 1;
shapfgdotshapfcRead = reshape(data(iInit:iEnd),[npf*npf,ngf,nd]);

iInit = iEnd + 1;
iEnd = iInit + npf*nfe - 1;
permRead = reshape(data(iInit:iEnd),[npf,nfe]);

iInit = iEnd + 1;
iEnd = iInit + paramLength - 1;
paramRead = reshape(data(iInit:iEnd),[paramLength,1]);

iInit = iEnd + 1;
iEnd = iInit + flagsLength - 1;
flagsRead = reshape(data(iInit:iEnd),[flagsLength,1]);

iInit = iEnd + 1;
iEnd = iInit + factorsLength - 1;
factorsRead = reshape(data(iInit:iEnd),[factorsLength,1]);

iInit = iEnd + 1;
iEnd = iInit + numBoundaries - 1;
bcmRead = reshape(data(iInit:iEnd),[numBoundaries,1]);

iInit = iEnd + 1;
iEnd = iInit + ncu*numBoundaries - 1;
bcsRead = reshape(data(iInit:iEnd),[ncu,numBoundaries]);

timeRead = data(iEnd+1);
