function [t1, t2] = cutgrid(p,t,a)

[~,c] = getcenter(p,t);
d = c(:,1)*a(1) + c(:,2)*a(2) + c(:,3)*a(3) + a(4); 
t1 = t(d<0,:);
t2 = t(d>0,:);


% figure(1); clf;
% hold on;
% nd = size(p,2);
% [ne,npf] = size(t1);    
% pn = zeros(npf,ne,nd);
% for k = 1:nd
%     pn(:,:,k) = reshape(p(t1,k),[ne npf])';
% end        
% patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),0*ones(npf,ne),'edgecolor','k');
% [ne,npf] = size(t2);    
% pn = zeros(npf,ne,nd);
% for k = 1:nd
%     pn(:,:,k) = reshape(p(t2,k),[ne npf])';
% end        
% patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),ones(npf,ne),'edgecolor','k');      
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])  
% 
