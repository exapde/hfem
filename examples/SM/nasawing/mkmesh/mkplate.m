function [p3, t3] = mkplate(p, t, th, cs)

[ne, npe] = size(t);
nd = size(p,2);

pn = zeros(npe,ne,nd);
for i = 1:nd
    pn(:,:,i) = reshape(p(t,i),[ne npe])';
end

nlg = normalvector(pn);
% averaging the normal vectors
[~,t1] = mkcgpt(permute(pn,[1 3 2]));
nsiz = max(t1(:));
n1   = zeros(nsiz,nd);
ndiv  = zeros(nsiz,1);
for i=1:ne
    il = t1(i,:);    
    n1(il,:) = n1(il,:) + reshape(nlg(:,i,:),[npe nd]);    
    ndiv(il) = ndiv(il) + 1;
end
ndiv = 1./ndiv;
n1 = bsxfun(@times,n1,ndiv);
nlg = zeros(npe,ne,nd);
for i = 1:nd
    nlg(:,:,i) = reshape(n1(t1,i),[ne npe])';
end

if cs==0
    pn1 = pn + (th/2)*nlg;
    pn2 = pn - (th/2)*nlg;
elseif cs==1
    pn1 = pn;
    pn2 = pn  + th*nlg;    
elseif cs==2
    pn1 = pn;
    pn2 = pn  - th*nlg; 
end


p3 = [reshape(pn1,[npe*ne nd]); reshape(pn2,[npe*ne nd])];
t1 = reshape(1:npe*ne,[npe ne])';
nsiz = npe*ne;
t3 = [t1 t1+nsiz];

[p3, t3] = fixmesh(p3, t3);

% nsiz = max(t(:));
% p1   = zeros(nsiz,nd);
% p2   = zeros(nsiz,nd);
% ndiv  = zeros(nsiz,1);
% for i=1:ne
%     il = t(i,:);    
%     p1(il,:) = p1(il,:) + reshape(pn1(:,i,:),[npe nd]);
%     p2(il,:) = p2(il,:) + reshape(pn2(:,i,:),[npe nd]);
%     ndiv(il) = ndiv(il) + 1;
% end
% ndiv = 1./ndiv;
% p1 = bsxfun(@times,p1,ndiv);
% p2 = bsxfun(@times,p2,ndiv);
% 
% p3 = [p1; p2];
% t3 = [t t+nsiz];

% nsiz = max(t(:));
% nl   = zeros(nsiz,nd);
% ndiv  = zeros(nsiz,1);
% for i=1:ne
%     il = t(i,:);    
%     nl(il,:) = nl(il,:) + reshape(nlg(:,i,:),[npe nd]);    
%     ndiv(il) = ndiv(il) + 1;
% end
% ndiv = 1./ndiv;
% nl = bsxfun(@times,nl,ndiv);


function nlg = normalvector(pn)

%nd = 3;
[npf,ne,nd] = size(pn);
porder = 1;
if npf==4
    elemtype=1;
    plocfc = [0  0; 1  0; 1  1; 0  1];
else
    elemtype=0;
    plocfc = [0  0; 1  0; 0  1];    
end
gpfc = plocfc;
npf = size(plocfc,1);
ngf = size(gpfc,1);

shapfc = mkshape(porder,plocfc,gpfc,elemtype);
dshapft  = reshape(permute(shapfc(:,:,2:nd),[2 3 1]),[ngf*(nd-1) npf]);
dpg = dshapft*reshape(pn,[npf ne*nd]);
dpg = permute(reshape(dpg,[ngf nd-1 ne nd]), [1 3 4 2]);    
dpg = reshape(dpg,[ngf*ne,nd,nd-1]);    

nlg = zeros(ngf*ne,nd);
nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
nlg   = bsxfun(@rdivide, nlg, jac);
nlg = reshape(nlg, [ngf ne nd]);


function [p,t]=fixmesh(p,t)
%FIXMESH  Remove duplicated/unused nodes and fix element orientation.
%   [P,T]=FIXMESH(P,T)

% Remove duplicated nodes:
snap=1e-8;
[foo,ix,jx]=unique(round(p/snap)*snap,'rows');
p=p(ix,:);
t=jx(t);
if size(t,2) == 1, t = t'; end  % This lines ensures the function works for one element

% Remove nodes that are not contained in t:
[pix,ix,jx]=unique(t);
t=reshape(jx,size(t));
p=p(pix,:);

