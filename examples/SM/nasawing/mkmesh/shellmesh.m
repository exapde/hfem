function [p3, t3] = shellmesh(p, t, thickness)

[ne, npe] = size(t);
nd = size(p,2);

pn = zeros(npe,ne,nd);
for i = 1:nd
    pn(:,:,i) = reshape(p(t,i),[ne npe])';
end

[pn1,pn2] = normalvector(pn, thickness');

nsiz = max(t(:));
p1   = zeros(nsiz,nd);
p2   = zeros(nsiz,nd);
ndiv  = zeros(nsiz,1);
for i=1:ne
    il = t(i,:);    
    p1(il,:) = p1(il,:) + reshape(pn1(:,i,:),[npe nd]);
    p2(il,:) = p2(il,:) + reshape(pn2(:,i,:),[npe nd]);
    ndiv(il) = ndiv(il) + 1;
end
ndiv = 1./ndiv;
p1 = bsxfun(@times,p1,ndiv);
p2 = bsxfun(@times,p2,ndiv);

p3 = [p1; p2];
t3 = [t t+nsiz];

% figure(1); clf;
% patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thickness','edgecolor','k')
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])
% 
% for i = 1:nd
%     pn(:,:,i) = reshape(p1(t,i),[ne npe])';
% end
% figure(2); clf;
% patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thickness','edgecolor','k')
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])
% 
% for i = 1:nd
%     pn(:,:,i) = reshape(p2(t,i),[ne npe])';
% end
% figure(3); clf;
% patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thickness','edgecolor','k')
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])
% 
% 
% % for i = 1:nd
% %     pn(:,:,i) = reshape(p3(t3,i),[ne npe])';
% % end
% % figure(2); clf;
% % patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thickness','edgecolor','k')
% % axis equal
% % axis tight
% % colormap(jet);
% % colorbar
% % view([75 27])
% 
% 
