function part = partition(mesh,master,app)

ne     = mesh.ne;

% nk = 5;            % default 10 blocks
% nk = round(ne/nk);      
% nk = 1:nk:ne;
% nb = [nk(1:end-1)' nk(2:end)'-1];
% nb(end) = ne;

ngrsiz = 1;
ngr    = ceil(ne/ngrsiz);
ngrne  = round(ne/ngr);      
nk = 1:ngrne:ne;
nb = [nk(1:end); [nk(2:end)-1,ne]];
part.nb    = nb;

for jj = 1:size(nb,2)
    e1 = nb(1,jj);
    e2 = nb(2,jj); 
    part.boug(jj) = bouprepro(master,app,mesh.t2t(:,e1:e2));
end



function bou = bouprepro(master,app,t2t)

neb = size(t2t,2);
nd = master.nd;
ng1 = master.ng1;

[sid,elem] = find(t2t<0);
t2t = reshape(t2t,[(nd+1)*neb,1]);
bct = app.bcm(-t2t((nd+1)*(elem-1)+sid))';
bcv = app.bcs(-t2t((nd+1)*(elem-1)+sid),:);
[bct,ix] = sort(bct);
sid = sid(ix);
elem = elem(ix);
bcv = bcv(ix,:);
nbcs = max(bct);
nbp = ones(nbcs+1,2);
for ibcs=1:nbcs
    nbp(ibcs,2) = nbp(ibcs,1)+numel(find(bct==ibcs))-1;
    nbp(ibcs+1,1) = nbp(ibcs,2)+1;
end
nbp = nbp(1:nbcs,:);
nft = numel(sid);

loc = reshape(1:(nd+1)*ng1,[ng1,nd+1]);
list = loc(:,sid);
list = reshape(bsxfun(@plus,reshape((elem-1)*(nd+1)*ng1,1,[]),list),[ng1*nft 1]);

bou.nft = nft;
bou.nbcs = nbcs;
bou.bcv = bcv;
bou.nbp = nbp;
bou.sid = sid;
bou.elem = elem;
bou.list = list;