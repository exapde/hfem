function [pn1,pn2,nlg] = normalvector(pn, tn)

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

%pn1 = pn;
tnd = repmat(tn/2,[ngf 1]);
tnd = repmat(tnd, [1 1 nd]);
pn1 = pn + tnd.*nlg;
pn2 = pn - tnd.*nlg;

% figure(2); clf;
% patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),tn,'edgecolor','k');
% hold on;
% patch(pn1(:,:,1),pn1(:,:,2),pn1(:,:,3),tn,'edgecolor','k');
% patch(pn2(:,:,1),pn2(:,:,2),pn2(:,:,3),tn,'edgecolor','k');
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])


