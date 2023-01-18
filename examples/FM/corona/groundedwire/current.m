function Ig = current(master,mesh,UDG,ib,Etmax)

% physical parameters
K = 1.5e-4; % (m^2/Vs)
e0 = 8.854e-12; % (F/m)

% lengthscale parameters  
%a = 0.01; % radius of the wire  (m)
L0 = 1;   % reference length scale (m)  

%Ec = 3.1*(1 + 0.0308*sqrt(1/a))*1e6; % Peek field
%Etmax = 40*1e3; % maximum background electric (V/m)
E0 = Etmax; % reference electric field (V/m)

% nondimensional parameters
Etstar = Etmax/E0;

%porder = mesh.porder;
perm = mesh.perm;
ne = mesh.ne;
nd = mesh.nd;
[npf,nfe] = size(perm);
ngf = master.ngf;

elcon = reshape(mesh.elcon,[npf nfe ne]);
shapft = squeeze(master.shapft(:,:,1));
dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);

in = find(mesh.f(:,end)==ib);
if isempty(in)
    error('Boundary is invalid.');
end

ns = length(in); 
Ig = 0;
for j = 1:ns
    i = in(j);
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
    i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
    j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
    pn = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));      
    udgg = shapft*UDG(perm(j1,i1),:,fi(1));                                 
    dpg   = dshapft*pn;
    dpg   = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 2]);  
                
    jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);
    nlg   = [dpg(:,2),-dpg(:,1)];
    nlg   = bsxfun(@rdivide, nlg, jac);
    ung   = udgg(:,3).*nlg(:,1) + (udgg(:,5)+Etstar).*nlg(:,2);  
    
    rhog = (e0*E0/L0)*udgg(:,2);
    Eng  = E0*abs(ung);
    
    Ig = Ig + sum(K*(rhog.*Eng).*jac.*master.gwfc);    
%     dgnodes(:,:,j) = pg;
%     u(:,:,j) = [udgg ung];        
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
