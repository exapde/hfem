
tol = 1e-6;

p = round([X Y Z]/1e-8)*1e-8;
comp = 0;  t0 = CQUAD(ID==comp,:); th0 = thickness(ID==comp); % Upper surface
comp = 1;  t1 = CQUAD(ID==comp,:); th1 = thickness(ID==comp); % Lower surface
comp = 2;  t2 = CQUAD(ID==comp,:); th2 = thickness(ID==comp); % Spars
comp = 3;  t3 = CQUAD(ID==comp,:); th3 = thickness(ID==comp); % Ribs

% break into smaller shells
[upt, upn] = uppersurfaces(p, t0, th0); % Upper surface
[lot, lon] = lowersurfaces(p, t1, th1); % Lower surface
[spt, spn] = spars(p, t2, th2, t3);     % Spars
[rbt, rbn] = ribs(p, t3, th3);          % Ribs

ts = rbt;
tn = rbn;
n = length(ts);
for i = 1:length(spt)
    n = n + 1;
    ts{n} = spt{i};
    tn(n) = spn(i);
end
for i = 1:length(upt)
    n = n + 1;
    ts{n} = upt{i};
    tn(n) = upn(i);
end
for i = 1:length(lot)
    n = n + 1;
    ts{n} = lot{i};
    tn(n) = lon(i);
end
tn = 2*tn;

% plot all small shells and their thickness
npf = 4; nd = 3;
figure(1); clf;
hold on;
for m = 1:length(ts)
    ne = size(ts{m},1);    
    pn = zeros(npf,ne,nd);
    for k = 1:nd
        pn(:,:,k) = reshape(p(ts{m},k),[ne npf])';
    end        
    patch(pn(:,:,1),pn(:,:,2),pn(:,:,3),tn(m)*ones(npf,ne),'edgecolor','k');
    view([75 27])        
    pause(0.25);
end
axis equal
axis tight
colormap(jet);
colorbar
view([75 27])        

% Create plates from shells
N = length(ts);
pp = cell(N,1);
tp = cell(N,1);
for i = 1:N
    [pp{i}, tp{i}] = mkplate(p, ts{i}, tn(i), 0);
end

% ind = [53 97 145 193];
ib{1} = [1 53];
ib{2} = [2 53 54];
ib{3} = [3 54 55];
ib{4} = [4 55 56];
ib{5} = [8 12 56 57];
ib{6} = [13 57 58];
for i = 7:44
    ib{i} = ib{i-1}+1;
end
ib{45} = [52 96];
ib{46} = [1 97];
ib{47} = [2 97 98];
ib{48} = [3 98 99];
ib{49} = [4 99 100];
ib{50} = [5 100 101];
ib{51} = [9 101 102];
for i = 52:93
    ib{i} = ib{i-1}+1;
end
ib{94} = [52 144];
ib{95} = [5 6 9];
ib{96} = [6 7 10];
ib{97} = [7 8 11];
for i = 1:length(ib)
    [pb{i}, tb{i}, ~, pp] = mkbeam2(p, pp, tp, ts, tn, ib{i}, tol);
end

% Beams on the top surface :
ib{98} = [145 146];

figure(2); clf;
hold on;
ind = 1:144;
for i = 1:length(ind)
    boundaryplot(pp{ind(i)},tp{ind(i)},{'facecolor','g','edgecolor','k','Linew',1});   
    pause(0.25);
end
ind = 1:length(pb);
for i = 1:length(ind)   
    boundaryplot(pb{ind(i)},tb{ind(i)},{'facecolor','r','edgecolor','k','Linew',1}); 
    pause(0.25);
end

return;

