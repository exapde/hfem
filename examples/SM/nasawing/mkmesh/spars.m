function [ts, thn] = spars(p, t2, th2, t3)

[f2, ~, l2, ~] = interedgemesh(p,t2,p,t3,1e-6);
for i = 1:length(l2)
    a = l2{i};
    b = f2(a,1:2);
    c = unique(p(b(:),:),'rows');
    d(i,:) = mean(c,1);
end
[~,ind]=sort(d(:,2));
e(1,:) = mean(d(ind(1:2),:),1);
e(2,:) = mean(d(ind(end-1:end),:),1);
x1 = e(1,1); y1 = e(1,2);
x2 = e(2,1); y2 = e(2,2);
s = (y2-y1)*d(:,1)-(x2-x1)*d(:,2)+x2*y1-x1*y2;
d1 = d(s>0,:); [~,i1] = sort(d1(:,2)); d1 = d1(i1,:);
d2 = d(s<0,:); [~,i2] = sort(d2(:,2)); d2 = d2(i2,:);

%figure(2);clf;plot(d1(:,1),d1(:,2),'o',d2(:,1),d2(:,2),'s',e(:,1),e(:,2),'*-');

N = 32;
thn = zeros(N,1);
ts = cell(N,1);
q2p = zeros(length(t3,1)); % Connectivity Quad 2 Plates`

[~,c] = getcenter(p, t2);
s = (y2-y1)*c(:,1)-(x2-x1)*c(:,2)+x2*y1-x1*y2;

n = size(d1,1)-1;
for i = 1:n
    m = d1(i,2)<c(:,2) & c(:,2)<d1(i+1,2) & s>0;    
    ts{i} = t2(m,:);    
    tm = th2(m);
    thn(i) = tm(1);
    % Connectivities quads to plates
    q2p(m) = i;
end
k = size(d2,1)-1;
for i = 1:k
    m = d2(i,2)<c(:,2) & c(:,2)<d2(i+1,2) & s<0;
    ts{n+i} = t2(m,:);
    tm = th2(m);
    thn(n+i) = tm(1);
end

% Update Connectivities quads to plates

for i=1:length(ts)
    for j=1:length(ts{i},1)
        q2p(ts{i}(j)) = i;
    end
end

return;

% npf = 4; nd = 3;
% figure(1); clf;
% hold on;
% for m = 1:length(ts)
%     ne = size(ts{m},1);    
%     pn = zeros(npf,ne,nd);
%     for k = 1:nd
%         pn(:,:,k) = reshape(p(ts{m},k),[ne npf])';
%     end        
%     patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thn(m)*ones(npf,ne),'edgecolor','k');        
%     pause
% end
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])        
% 



