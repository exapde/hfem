
tol = 1e-6;

p = round([X Y Z]/1e-8)*1e-8;
comp = 2;  t2 = CQUAD(ID==comp,:); th2 = thickness(ID==comp); % Spars
comp = 3;  t3 = CQUAD(ID==comp,:); th3 = thickness(ID==comp); % Ribs

% Break the shell surface into shells
shell2 = breakshells(p, t2, th2, 0); % Spars
shell3 = breakshells(p, t3, th3, 0); % Ribs

N = 13;
thn = zeros(N,1);
ts = cell(N,1);
ts2 = shell2{3};
ts3{1} = shell3{1}{1};
for i=2:1:(length(shell3{1})-1)
    ts3{i}=shell3{1}{i+1};
end
ind = [1 2 3 6 5 4 7 8];
for i=1:8
    ts{i} = t2(ts2{ind(i)},:);
    thn(i) = th2(ts2{ind(i)}(1));
end
ind = [1 2 3 4 5];
for i=1:5
    ts{8+i} = t3(ts3{ind(i)},:);
    thn(8+i) = th3(ts3{ind(i)}(1));
end
thn = 10*thn;

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
end
axis equal
axis tight
colorbar('FontSize',14);
set(gca,'FontSize',16);
colormap(jet);
view([40 25])    

% Create plates from shells
pp = cell(N,1);
tp = cell(N,1);
cp = [1 1 1 1 2 2 2 2 0 0 0 0 0];
for i = 1:N
    [pp{i}, tp{i}] = mkplate(p, ts{i}, thn(i), cp(i));
end

% Create beams from shells and fix plates from beams
cb = [2 0 0 0 2 3 1 1 1 3];
sb{1} = [1 9]; 
sb{2} = [1 2 10];
sb{3} = [2 3 11];
sb{4} = [3 4 12];
sb{5} = [4 13]; 
sb{6} = [5 13]; 
sb{7} = [5 6 12]; 
sb{8} = [6 7 11]; 
sb{9} = [7 8 10]; 
sb{10} = [8 9]; 
for i = 1:length(cb)
    [pb{i}, tb{i}, eb{i}, pp] = mkbeamplate(p, ts, pp, thn, sb{i}, cb(i), tol);
end

figure(2); clf;
hold on;
ind = 1:length(cb);
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},{'facecolor','r','edgecolor','k','Linew',1}); 
end
ind = 1:length(cp);
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},{'facecolor','g','edgecolor','k','Linew',1}); 
end
axis equal
axis tight
set(gca,'FontSize',16);
colormap(jet);
view([40 25]);    



