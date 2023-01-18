function [dgnodes,u] = normalelectricfield(master,mesh,UDG,UH,ib,tau)

%porder = mesh.porder;
perm = mesh.perm;
ne = mesh.ne;
nd = mesh.nd;
[npf,nfe] = size(perm);

elcon = reshape(mesh.elcon,[npf nfe ne]);
% shapft = squeeze(master.shapft(:,:,1));
% ngf = master.ngf;
% dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
xi     = linspace(0,1,40)';
shapmf = mkshape(master.porder,master.plocfc,xi,1);
shapft = shapmf(:,:,1)';
ngf = size(shapft,1);
dshapft  = reshape(permute(shapmf(:,:,2:nd),[2 3 1]),[ngf*(nd-1) npf]);

in = find(mesh.f(:,end)==ib);
if isempty(in)
    error('Boundary is invalid.');
end

nch = size(UH,1);
nc = size(UDG,2);
nf = size(mesh.f,1);
UH = permute(reshape(UH,[nch,npf,nf]),[2 1 3]);

ns = length(in); 
dgnodes = zeros(ngf,nd,ns);
u = zeros(ngf,nc+1,ns);
for j = 1:ns
    i = in(j);
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
    i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
    j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
    pn = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));  
    pg = shapft*pn;
    udgg = shapft*UDG(perm(j1,i1),:,fi(1)); 
    uhg = shapft*UH(:,:,i); 
    dpg   = dshapft*pn;
    dpg   = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 2]);  
                
    if nd==2
        jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
        ung   = udgg(:,3).*nlg(:,1) + (udgg(:,5)+1).*nlg(:,2) + tau*(udgg(:,2)-uhg(:,2));  
    elseif nd==3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
        ung   = udgg(:,1).*nlg(:,1) + udgg(:,2).*nlg(:,2) + udgg(:,3).*nlg(:,3);  
    end                     
    
    dgnodes(:,:,j) = pg;
    u(:,:,j) = [udgg ung];        
end


% figure(1); clf;
% hold on;
% for k = 1:ns
%     plot3(dgnodes(:,1,k),dgnodes(:,2,k),u(:,1,k),'-k');
% end
% hold off;
% axis normal;
% axis tight;
% view(3);
% 
% figure(2); clf;
% hold on;
% for k = 1:ns
%     plot3(dgnodes(:,1,k),dgnodes(:,2,k),u(:,2,k),'-k');
% end
% hold off;
% axis normal;
% axis tight;
% view(3);
% 
% figure(3); clf;
% hold on;
% for k = 1:ns
%     plot3(dgnodes(:,1,k),dgnodes(:,2,k),sqrt(u(:,1,k).^2+u(:,2,k).^2),'-k');
% end
% hold off;
% axis normal;
% axis tight;
% view(3);
