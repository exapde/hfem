function [w,fj] = getfbou(mesh,master,E,E0,Ecstar,delta_peek)

% Et = Etstar;     
% E(:,1,:) = UDG(:,3,:); E(:,2,:) = UDG(:,5,:)+Et;
[dgnodes,u] = normalelectricfield(master,mesh,E,-5,[1 2],(1-delta_peek)*Ecstar);
x = dgnodes(:,1,:);
y = dgnodes(:,2,:);
f = u(:,3,:);
q = u(:,4,:);
% rho = f-((1-delta_peek)*Ecstar);
% ind = find(f<(1-delta_peek)*Ecstar);
% rho(ind) = 0; 
x=x(:); y=y(:); f=f(:); q = q(:);
t = cart2pol(x,y);
i = find(t<0);
t(i) = 2*pi+t(i);
[tj,jj] = sort(t);
w = trapz(tj,q(jj));

[~,ii] = sort(abs(f-(1-delta_peek)*Ecstar));
[~,i1]=sort(t(ii(1:8)));
ii = ii(i1([1 end]));
fj = f(jj)*E0/1e6;
qj = q(jj)*E0/1e6;
t1 = t(ii(1));
t2 = t(ii(2));

figure(1);clf;
fill([0 2*pi 2*pi 0],[(1-delta_peek)*E0*Ecstar/1e6 (1-delta_peek)*E0*Ecstar/1e6 (1+delta_peek)*E0*Ecstar/1e6 (1+delta_peek)*E0*Ecstar/1e6],[0.4 0.4 0.4]*1.5);
hold on;
fill([t1 t2 t2 t1],[min(fj) min(fj) max(fj) max(fj)],[0.4 0.4 0.4]*1.5);
plot(tj,fj,'b-',tj,E0*Ecstar*ones(size(t))/1e6,'--k','LineWidth',1.5);
plot(tj,(1+delta_peek)*E0*Ecstar*ones(size(t))/1e6,'-k','LineWidth',1.5);
plot(tj,(1-delta_peek)*E0*Ecstar*ones(size(t))/1e6,'-k','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));

figure(2);clf;
fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(qj) min(qj) max(qj) max(qj)],[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,qj,'-b','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));

figure(3);clf;
hold on;
plot(tj,qj*1e6/E0,'-b','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));


function [dgnodes,u,w] = normalelectricfield(master,mesh,UDG,ib,ind,f0)

%porder = mesh.porder;
perm = mesh.perm;
ne = mesh.ne;
nd = mesh.nd;
[npf,nfe] = size(perm);

elcon = reshape(mesh.elcon,[npf nfe ne]);
shapft = squeeze(master.shapft(:,:,1));
ngf = master.ngf;
dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
% xi     = linspace(0,1,40)';
% shapmf = mkshape(master.porder,master.plocfc,xi,1);
% shapft = shapmf(:,:,1)';
% ngf = size(shapft,1);
% dshapft  = reshape(permute(shapmf(:,:,2:nd),[2 3 1]),[ngf*(nd-1) npf]);

in = find(mesh.f(:,end)==ib);
if isempty(in)
    error('Boundary is invalid.');
end

w = 0;
ns = length(in); 
dgnodes = zeros(ngf,nd,ns);
u = zeros(ngf,length(ind)+2,ns);
for j = 1:ns
    i = in(j);
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
    i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
    j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
    pn = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));  
    pg = shapft*pn;
    udgg = shapft*UDG(perm(j1,i1),ind,fi(1));                                 
    dpg   = dshapft*pn;
    dpg   = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 2]);  
                
    if nd==2
        jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
        ung   = udgg(:,1).*nlg(:,1) + udgg(:,2).*nlg(:,2);  
    elseif nd==3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
        ung   = udgg(:,1).*nlg(:,1) + udgg(:,2).*nlg(:,2) + udgg(:,3).*nlg(:,3);  
    end                     

    f = abs(ung);
    q = f-f0;    
    q(f<f0) = 0;    
    dgnodes(:,:,j) = pg;
    u(:,:,j) = [udgg f q];        
    
    w = w + sum(q.*jac.*master.gwfc);
end

