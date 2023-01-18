function [ts, thn] = ribs(p, t3, th3)

% Break the shell surface into shells
shell3 = breakshells(p, t3, th3, 1); % Ribs

% Get center of shells
center3 = cell(length(shell3),1);
for i=1:length(shell3)
    for j = 1:length(shell3{i})
        center3{i}(j,:) = getcenter(p, t3(shell3{i}{j},:));
    end
end

N = 32;
thn = zeros(N,1);
ts = cell(N,1);
n = 0;

[~,ind] = sort(center3{1}(:,2));
for j=1:length(ind)
    n = n + 1;
    i = ind(j);    
    ts{n} = t3(shell3{1}{i},:);
    thn(n) = th3(shell3{1}{i}(1));    
end

[t1, t2] = cutgrid(p,ts{23},[0 1 0 -21.8]);
thn23 = thn(23);
ts{23} = ts{22}; thn(23) = thn(22);
ts{22} = ts{21}; thn(22) = thn(21); 
ts{21} = ts{20}; thn(21) = thn(20);
ts{20} = t1;     thn(20) = thn23;
ts{49} = t2;     thn(49) = thn23;

[t1, t2] = cutgrid(p,ts{36},[0 1 0 -21.8]);
for i = 48:-1:39
    ts{i} = ts{i-1};
    thn(i) = thn(i-1);
end
ts{36} = t1;
ts{38} = t2;
thn(38)  = thn(36);

% break ts{5}
for j = 6:9
    [f2, ~, l2, ~] = interedgemesh(p,ts{5},p,ts{j},1e-6);    
    a = l2{1};
    b = f2(a,1:2);
    c = unique(p(b(:),:),'rows');
    d(j-5,:) = mean(c,1);    
end
[~,c] = getcenter(p, ts{5});
m = c(:,1)<d(1,1);  
tq{1} = ts{5}(m,:);        
thq(1) = thn(5);
for i = 1:3
    m = d(i,1)<c(:,1) & c(:,1)<d(i+1,1);  
    tq{i+1} = ts{5}(m,:);        
    thq(i+1) = thn(5);
end

for i = 52:-1:9
    ts{i} = ts{i-3};
    thn(i) = thn(i-3);
end
for i = 5:8
    ts{i} = tq{i-4};
    thn(i) = thq(i-4);
end

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
%     view([75 27])            
% end
% axis equal
% axis tight
% colormap(jet);
% colorbar
% view([75 27])        
