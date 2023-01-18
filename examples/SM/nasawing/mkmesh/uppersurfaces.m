function [ts, thn, thm, q2p] = uppersurfaces(p, t3, th3)

% Break the shell surface into shells
shell3 = breakshells(p, t3, th3, 1);

% Get center of shell panels
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

% Create Connectivities and thickness for each shell panels.
[~,ind] = sort(center3{1}(:,2));
nquads = length(ind);
for j=1:nquads
    n = n + 1;
    i = ind(j);    
    ts{n} = t3(shell3{1}{i},:);
    thn(n) = th3(shell3{1}{i}(1));        
end

[t1, t2] = cutgrid(p,ts{15},[0 1 0 -9.5]);
thn15 = thn(15);
ii = [11 20];
for i = 15:-1:(ii(1)+1)
    ts{i} = ts{i-1};
    thn(i) = thn(i-1);
end
ts{ii(1)} = t1;
thn(ii(1)) = thn15;
for i = (n+1):-1:(ii(2)+1)
    ts{i} = ts{i-1};
    thn(i) = thn(i-1);
end
ts{ii(2)} = t2;
thn(ii(2)) = thn15;

for m = 1:length(ts)
    thm(m) = size(ts{m},1);
end

% Update Connectivities quads to plates
q2p = zeros(length(t3,1));
for i=1:length(shell3{1})
    for j=1:length(shell3{1}{i})
        q2p(shell3{1}{i}(j)) = i;
    end
end

return;




npf = 4; nd = 3;
figure(1); clf;
hold on;
for m = 1:length(ts)
    ne = size(ts{m},1);    
    pn = zeros(npf,ne,nd);
    for k = 1:nd
        pn(:,:,k) = reshape(p(ts{m},k),[ne npf])';
    end        
    patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thn(m)*ones(npf,ne),'edgecolor','k');
    view([75 27])        
    pause
end
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])        
