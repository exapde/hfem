function [dgnodes,u] = potentialfield(master,mesh,UDG,ib,ind)

%porder = mesh.porder;
perm = mesh.perm;
ne = mesh.ne;
nd = mesh.nd;
[npf,nfe] = size(perm);

elcon = reshape(mesh.elcon,[npf nfe ne]);

%shapft = squeeze(master.shapft(:,:,1));
xi     = linspace(0,1,40)';
shapmf = mkshape(master.porder,master.plocfc,xi,1);
shapft = shapmf(:,:,1)';
ngf = size(shapft,1);

in = find(mesh.f(:,end)==ib);
if isempty(in)
    error('Boundary is invalid.');
end

ns = length(in); 
dgnodes = zeros(ngf,nd,ns);
u = zeros(ngf,length(ind),ns);
for j = 1:ns
    i = in(j);
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
    i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
    j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
    pn = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));  
    pg = shapft*pn;
    udgg = shapft*UDG(perm(j1,i1),ind,fi(1));                                     
    dgnodes(:,:,j) = pg;
    u(:,:,j) = udgg;        
end

% figure(1); clf;
% hold on;
% for k = 1:ns
%     plot(dgnodes(:,2,k),2.5e4*u(:,1,k),'-k');
% end
% hold off;
% axis normal;
% axis tight;
%view(3);

% figure(2); clf;
% hold on;
% for k = 1:ns
%     plot3(dgnodes(:,1,k),dgnodes(:,2,k),u(:,1,k),'-k');
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
