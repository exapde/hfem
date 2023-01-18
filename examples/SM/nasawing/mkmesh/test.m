

%scl = [25 25 50 50];
scl = [1 1 1 1];
for comp = 0:3
scale = scl(comp+1);
[p, t] = shellmesh(round([X Y Z]/1e-8)*1e-8, CQUAD(ID==comp,:), scale*thickness(ID==comp));

face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]]';
[npf,nfe] = size(face);
[ne,npe] = size(t);
nd = size(p,2);

p1 = zeros(npf,nfe,ne,nd);
for i = 1:nd
    p1(:,:,:,i) = permute(reshape(p(t(:,face(:)),i),[ne npf nfe]),[2 3 1]);
end

tn = repmat(scale*thickness(ID==comp)',[nfe 1]);
figure(comp+1); clf;
patch(reshape(p1(:,:,:,1),npf,[]),reshape(p1(:,:,:,2),npf,[]),reshape(p1(:,:,:,3),npf,[]),tn(:)','edgecolor','k')
axis equal
axis tight
colormap(jet);
% colorbar
view([75 27]);

fn = ['comp' num2str(comp+4)];
print('-dpng',fn);
end

p = round([X Y Z]/1e-8)*1e-8;
t1 = CQUAD(ID==2,:);
t2 = CQUAD(ID==3,:);
[f1, f2] = interedgemesh(p,t1,p,t2,1e-8);


% tn = repmat(scale*thickness(ID==2)',[nfe 1]);
% figure(2); clf;
% patch(reshape(p1(:,:,ID==2,1),npf,[]),reshape(p1(:,:,ID==2,2),npf,[]),reshape(p1(:,:,ID==2,3),npf,[]),tn(:)','edgecolor','k')
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])


% % Allocate nodes
% dgnodes=zeros(npl,nd,nt);
% for dim=1:nd
%   for node=1:npv
%     dp=philocal(:,node)*p(t(:,node),dim)';
%     dgnodes(:,dim,:)=dgnodes(:,dim,:)+permute(dp,[1,3,2]);
%   end
% end


% porder = 1;
% elemtype = 1;
% nodetype = 0;
% bndexpr={'true'};
% mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

% nd = 3;
% elemtype=1;
% porder = 1;
% plocfc = [0  0; 1  0; 1  1; 0  1];
% % gpfc =  [0.211324865405187   0.211324865405187;
% %          0.788675134594813   0.211324865405187;         
% %          0.788675134594813   0.788675134594813;
% %          0.211324865405187   0.788675134594813; 
% %          0.5 0.5];   
% gpfc = plocfc;
% npf = size(plocfc,1);
% ngf = size(gpfc,1);
% ne = size(CQUAD(ID==3,:),1);
% 
% pn = zeros(npf,ne,nd);
% pn(:,:,1)=X(CQUAD(ID==3,:))';
% pn(:,:,2)=Y(CQUAD(ID==3,:))';
% pn(:,:,3)=Z(CQUAD(ID==3,:))';
% tn=thickness(ID==3)';
% 
% snap = 1e-8;
% pn = round(pn/snap)*snap;
% tn = round(tn/snap)*snap;
% 
% shapfc = mkshape(porder,plocfc,gpfc,elemtype);
% dshapft  = reshape(permute(shapfc(:,:,2:nd),[2 3 1]),[ngf*(nd-1) npf]);
% dpg = dshapft*reshape(pn,[npf ne*nd]);
% dpg = permute(reshape(dpg,[ngf nd-1 ne nd]), [1 3 4 2]);    
% dpg = reshape(dpg,[ngf*ne,nd,nd-1]);    
% 
% nlg = zeros(ngf*ne,nd);
% nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
% nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
% nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
% jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
% nlg   = bsxfun(@rdivide, nlg, jac);
% 
% nlg = permute(reshape(nlg, [ngf ne nd]), [1 3 2]);
% 
