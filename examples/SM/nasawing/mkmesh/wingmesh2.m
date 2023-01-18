
tol = 1e-6;

p = round([X Y Z]/1e-8)*1e-8;
comp = 0;  t0 = CQUAD(ID==comp,:); th0 = thickness(ID==comp); % Upper surface
comp = 1;  t1 = CQUAD(ID==comp,:); th1 = thickness(ID==comp); % Lower surface
comp = 2;  t2 = CQUAD(ID==comp,:); th2 = thickness(ID==comp); % Spars
comp = 3;  t3 = CQUAD(ID==comp,:); th3 = thickness(ID==comp); % Ribs

% Break the shell surface into shells
shell0 = breakshells(p, t0, th0, 1); % Upper surface
shell1 = breakshells(p, t1, th1, 1); % Lower surface
shell2 = breakshells(p, t2, th2, 0); % Spars
shell3 = breakshells(p, t3, th3, 0); % Ribs

npf = 4; nd=3; ne = length(shell1{1}{34});
pn = zeros(npf,ne,nd);
for i = 1:nd
    pn(:,:,i) = reshape(p(t1(shell1{1}{34},:),i),[ne npf])';
end
pn = reshape(mean(pn,1),[ne nd]);
shell1{1}{34} = shell1{1}{34}(pn(:,2)<5);

N = 21;
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
ind = [33 40 42 45];
for i=1:length(ind)
    ts{13+i} = t0(shell0{1}{ind(i)},:);
    thn(13+i) = th0(shell0{1}{ind(i)}(1));
end
ind = [27 29 34 42]; % 34 42];
for i=1:length(ind)
    ts{17+i} = t1(shell1{1}{ind(i)},:);
    thn(17+i) = th1(shell1{1}{ind(i)}(1));
end

% make the thickness 10 times larger
thn = 1*thn;  % make the thickness 10 times larger

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
colormap(jet);
colorbar
view([75 27])        

% Create plates from shells
pp = cell(N,1);
tp = cell(N,1);
cp = [1 1 1 1 2 2 2 2 0 0 0 0 0 1 1 1 1 2 2 2 2];
for i = 1:N
    [pp{i}, tp{i}] = mkplate(p, ts{i}, thn(i), cp(i));
end

% Create beams from shells and fix plates from beams
cb = [2 0 0 0 2 3 1 1 1 3 2 0 0 0 2 3 1 1 1 3];
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
sb{11} = [14 9]; 
sb{12} = [14 15 10]; 
sb{13} = [15 16 11];
sb{14} = [16 17 12];
sb{15} = [17 13]; 
sb{16} = [21 13]; 
sb{17} = [21 20 12]; 
sb{18} = [20 19 11];
sb{19} = [19 18 10];
sb{20} = [18 9]; 
for i = 1:length(cb)
    [pb{i}, tb{i}, eb{i}, pp] = mkbeamplate(p, ts, pp, thn, sb{i}, cb(i), tol);
end

% Create beams from two plates
I = [1 14; 2 15; 3 16; 4 17; 1 18; 2 19; 3 20; 4 21; 8 14; 7 15; 6 16; 5 17; 8 18; 7 19; 6 20; 5 21];
for k = 1:size(I,1)
    i = I(k,1); j = I(k,2);
    [pb{20+k},tb{20+k}] = mkbeamfrom2plates(pp{i}, tp{i}, thn(i), pp{j}, tp{j}, thn(j));
end

% Create cells from beams
sc{1} = [1 11 21];
sc{2} = [2 12 21 22];
sc{3} = [3 13 22 23];
sc{4} = [4 14 23 24];
sc{5} = [5 15 24];
sc{6} = [1 20 25];
sc{7} = [2 19 25 26];
sc{8} = [3 18 26 27];
sc{9} = [4 17 27 28];
sc{10} = [5 16 28];
sc{11} = [10 11 29];
sc{12} = [9 12 29 30];
sc{13} = [8 13 30 31];
sc{14} = [7 14 31 32];
sc{15} = [6 15 32];
sc{16} = [10 20 33];
sc{17} = [9 19 33 34];
sc{18} = [8 18 34 35];
sc{19} = [7 17 35 36];
sc{20} = [6 16 36];
for i = 1:length(sc)
    [pc{i},tc{i}] = mkcellfrombeams(pb,tb,sc{i});
end

figure(2); clf;
hold on;
ind = 1:length(tb);
for i = 1:length(ind)
    boundaryplot(pb{ind(i)},tb{ind(i)},{'facecolor','r','edgecolor','k','Linew',1}); 
end
ind = 1:length(tp);
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},{'facecolor','g','edgecolor','k','Linew',1}); 
end
ind = 1:length(tc);
for i = 1:length(ind)
    boundaryplot(pc{ind(i)},tc{ind(i)},{'facecolor','b','edgecolor','k','Linew',1}); 
end
axis equal
axis tight
set(gca,'FontSize',16);
colormap(jet);
view([40 25])    

% Connect plates, beams, and cells to make the 3D wing mesh 
[pw,tw] = connectmesh(pp{1},tp{1},pp{2},tp{2});
for i = 3:length(tp)
    [pw,tw] = connectmesh(pw,tw,pp{i},tp{i});
end
for i = 1:length(tb)
    [pw,tw] = connectmesh(pw,tw,pb{i},tb{i});
end
for i = 1:length(tc)
    [pw,tw] = connectmesh(pw,tw,pc{i},tc{i});
end
figure(3); clf;
boundaryplot(pw,tw,{'facecolor','g','edgecolor','k','Linew',1}); 
axis equal
axis tight
set(gca,'FontSize',16);
colormap(jet);
view([40 25])    

% figure(4); clf;
% hold on;
% for m = 1:5
%     ne = size(ts{m},1);    
%     pn = zeros(npf,ne,nd);
%     for k = 1:nd
%         pn(:,:,k) = reshape(p(ts{m},k),[ne npf])';
%     end        
%     patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),thn(m)*ones(npf,ne),'edgecolor','k');    
% end
% for m = 14:17
%     ne = size(ts{m},1);    
%     pn = zeros(npf,ne,nd);
%     for k = 1:nd
%         pn(:,:,k) = reshape(p(ts{m},k),[ne npf])';
%     end        
%     patch(pn(:,:,1),pn(:,:,2),pn(:,:,3)+1,thn(m)*ones(npf,ne),'edgecolor','k');    
% end
% for m = 18:21
%     ne = size(ts{m},1);    
%     pn = zeros(npf,ne,nd);
%     for k = 1:nd
%         pn(:,:,k) = reshape(p(ts{m},k),[ne npf])';
%     end        
%     patch(pn(:,:,1),pn(:,:,2),pn(:,:,3)-1,thn(m)*ones(npf,ne),'edgecolor','k');    
% end
% axis equal
% axis tight
% set(gca,'FontSize',16);
% colormap(jet);
% view([40 25])    
