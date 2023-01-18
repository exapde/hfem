function [shells, ns, es] = breakshells(p, t, thickness, cs)

if nargin<4
    cs = 0;
end

if cs==0
[ne, npf] = size(t);
nd = size(p,2);

pn = zeros(npf,ne,nd);
for i = 1:nd
    pn(:,:,i) = reshape(p(t,i),[ne npf])';
end
[~,~,nl] = normalvector(pn, 0.01*ones(1,ne));
nl = reshape(mean(nl,1),[ne nd]);

tol = 2e-2;
nl = round(nl/tol)*tol;
[ns,~,ic] = unique(nl,'rows');    
j = max(ic);
for i = 1:j
    es{i} = find(ic==i);
    tn = thickness(es{i});
    tm = unique(tn);
    k = length(tm);
    for m = 1:k                
        shells{i}{m} = es{i}(tn==tm(m));
    end        
end
else
    ns = [];
    es = [];
    tm = unique(thickness);
    k = length(tm);
    for m = 1:k                
        shells{1}{m} = find(thickness==tm(m));
    end            
end

% [nl,ind] = sortrows(nl);
% tol = 5.0e-2;
% j = 1;
% ns = nl;
% es = cell(ne,1);
% es{j} = 1;
% for i=2:ne
%     if norm(nl(i,:)-ns(j,:))>tol        
%         j = j + 1;
%         ns(j,:) = nl(i,:);
%         es{j} = i;
%     else
%         es{j} = [es{j} i];
%     end
% end
% ns = ns(1:j,:);
% nes = zeros(1,j);
% for i = 1:j
%     es{i} = ind(es{i});
%     nes(i) = length(es{i});
%     ns(i,:) = mean(nl(es{i},:),1);
% end
% 
% tol = 0.3;
% for i=1:j
%     ne = length(es{i});    
%     ts = t(es{i},:);
%     pn = zeros(npf,ne,nd);
%     for k = 1:nd
%         pn(:,:,k) = reshape(p(ts,k),[ne npf])';
%     end
%     pn = reshape(mean(pn,1),[ne nd]);
%     d = pn*(ns(i,:)');   
%     d = round(d/tol)*tol;    
%     [~,~,ic] = unique(d);    
%     k = max(ic);    
%     for m = 1:k
%         shells{i}{m} = es{i}(ic==m);
%     end
% end

isplot = 0;
if isplot==1
[ne, npf] = size(t);
pn = zeros(npf,ne,nd);
for i = 1:nd
    pn(:,:,i) = reshape(p(t,i),[ne npf])';
end
figure(1); clf;
patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thickness','edgecolor','k')
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])

% for i=1:j
%     ne = length(es{i});    
%     ts = t(es{i},:);
%     pn = zeros(npf,ne,nd);
%     for k = 1:nd
%         pn(:,:,k) = reshape(p(ts,k),[ne npf])';
%     end
%     %size(thickness(es{i-1})')
%     figure(i+1); clf;
%     patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thickness(es{i})','edgecolor','k')
%     axis equal
%     axis tight
%     colormap(jet);
%     colorbar
%     view([75 27])    
% end

for i=1:j
    figure(i+1); clf;
    hold on;
    for m = 1:length(shells{i})
        ne = length(shells{i}{m});    
        ts = t(shells{i}{m},:);
        pn = zeros(npf,ne,nd);
        for k = 1:nd
            pn(:,:,k) = reshape(p(ts,k),[ne npf])';
        end        
        patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thickness(shells{i}{m})','edgecolor','k');
    end
    axis equal
    axis tight
    colormap(jet);
    colorbar
    view([75 27])        
end

% figure(1); clf;
% patch('vertices',p,'faces',t,'FaceVertexCData',thickness,'edgec','none');
% hold on;
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])

end

return;




comp=2; [shell2, ns2, es2] = breakshells(round([X Y Z]/1e-8)*1e-8, CQUAD(ID==comp,:), thickness(ID==comp));
comp=3; [shell3, ns3, es3] = breakshells(round([X Y Z]/1e-8)*1e-8, CQUAD(ID==comp,:), thickness(ID==comp));

ts2 = shell2{3};
ts3{1} = shell3{1}{1};
for i=2:1:(length(shell3{1})-1)
    ts3{i}=shell3{1}{i+1};
end

p = round([X Y Z]/1e-8)*1e-8;
npf = 4; nd = 3;
comp = 2;  t = CQUAD(ID==comp,:); tn = thickness(ID==comp);
figure(1); clf;
hold on;
for m = 1:length(ts2)
    ne = length(ts2{m});    
    ts = t(ts2{m},:);
    pn = zeros(npf,ne,nd);
    for k = 1:nd
        pn(:,:,k) = reshape(p(ts,k),[ne npf])';
    end        
    patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),tn(ts2{m})','edgecolor','k');
end
comp = 3;  t = CQUAD(ID==comp,:); tn = thickness(ID==comp);
for m = 1:length(ts3)
    ne = length(ts3{m});    
    ts = t(ts3{m},:);
    pn = zeros(npf,ne,nd);
    for k = 1:nd
        pn(:,:,k) = reshape(p(ts,k),[ne npf])';
    end        
    patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),tn(ts3{m})','edgecolor','k');
end
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])        

scale = 20;
face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]]';
[npf,nfe] = size(face); nd = 3;
comp = 2;  t = CQUAD(ID==comp,:); tn = thickness(ID==comp);
figure(1); clf;
hold on;
for m = 1:length(ts2)    
    [pm, tm] = shellmesh(p, t(ts2{m},:), scale*tn(ts2{m}));
    
    [ne,npe] = size(tm);
    p1 = zeros(npf,nfe,ne,nd);
    for i = 1:nd
        p1(:,:,:,i) = permute(reshape(pm(tm(:,face(:)),i),[ne npf nfe]),[2 3 1]);
    end
    
    th = repmat(scale*tn(ts2{m})',[nfe 1]);    
    patch(reshape(p1(:,:,:,1),npf,[]),reshape(p1(:,:,:,2),npf,[]),reshape(p1(:,:,:,3),npf,[]),th(:)','edgecolor','k')
end
comp = 3;  t = CQUAD(ID==comp,:); tn = thickness(ID==comp);
for m = 1:length(ts3)    
    [pm, tm] = shellmesh(p, t(ts3{m},:), scale*tn(ts3{m}));
    
    [ne,npe] = size(tm);
    p1 = zeros(npf,nfe,ne,nd);
    for i = 1:nd
        p1(:,:,:,i) = permute(reshape(pm(tm(:,face(:)),i),[ne npf nfe]),[2 3 1]);
    end
    
    th = repmat(scale*tn(ts3{m})',[nfe 1]);    
    patch(reshape(p1(:,:,:,1),npf,[]),reshape(p1(:,:,:,2),npf,[]),reshape(p1(:,:,:,3),npf,[]),th(:)','edgecolor','k')
end
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])        


comp = 2;  t2 = CQUAD(ID==comp,:); tn2 = thickness(ID==comp);
comp = 3;  t3 = CQUAD(ID==comp,:); tn3 = thickness(ID==comp);

ta2 = ts2{1};
for i = 2:length(ts2)
    ta2 = [ta2; ts2{i}];
end
ta3 = ts3{1};
for i = 2:length(ts3)
    ta3 = [ta3; ts3{i}];
end

% [f1, f2, l1, l2] = interedgemesh(p,t2(ta2,:),p,t3(ta3,:),1e-6);
% [f1, f2, n1, n2] = interedgemesh(p,t2(ts2{1},:),p,t3(ts3{1},:),1e-6);

[f1, f2, n1, n2, c1, c2] = intersectshells(p,t2(ts2{1},:),p,t3(ts3{1},:),1e-6);

beamline = p([f1(:,1); f1(end,2)],:);
normal = n1;
normal(:,:,2) = n2;
thickness = 50*[tn2(ts2{1}(1)) tn3(ts3{1}(1))];
shellcenter = [c1; c2];
[pa, ta, ea] = mkbeam(beamline, normal, thickness, shellcenter, 2);
pars={'facecolor','r','edgecolor','k','Linew',1};
boundaryplot(pa,ta,pars); view(3);

plot3(pa(:,1),pa(:,2),pa(:,3),'ro','LineWidth',1.5);

face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]]';
[npf,nfe] = size(face);
[ne,npe] = size(ta);
nd = size(pa,2);
p1 = zeros(npf,nfe,ne,nd);
for i = 1:nd
    p1(:,:,:,i) = permute(reshape(pa(ta(:,face(:)),i),[ne npf nfe]),[2 3 1]);
end
c = zeros(nfe*ne,1);
figure(2); clf;
patch(reshape(p1(:,:,:,1),npf,[]),reshape(p1(:,:,:,2),npf,[]),reshape(p1(:,:,:,3),npf,[]),c,'edgecolor','k')
axis equal
axis tight
colormap(jet);
% colorbar
view([75 27]);

pars={'facecolor','g','edgecolor','k','Linew',1};
figure(3);clf; boundaryplot(pa,ta,pars); view(3);



